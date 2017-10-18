import math
import memfiles
import parseopt
import random
import sets
import sequtils
import strutils
import tables
import terminal

const VERSION = "0.0.2, 12 Oct 2017"
const AUTHOR = "Mark Pinese <m.pinese@garvan.org.au>"


# TODO: Add warning somewhere for AF mismatch
# TODO: Special rules for resampling:
#   * Sampling should stay within xsome "class" (classes: X, Y, MT, autosomes)
#   * Two choices of AF: dosages (empirical) or model
#   * Special treatment of missing loci (in model but not in dosages).  
#     These should *not* be resampled, but rather removed and added to offset.
#     One hacky way to achieve this is to change the missing locus ID to an
#     ID which is definitely not in the dosages data.
#   * Try and do it without replacement


type
  Dosages = tuple[
    samples: seq[string],
    vids: seq[string],
    vid2idx: Table[string, int],
    afs: seq[float],
    afbins: seq[int],
    afbin2idx: seq[seq[int]],
    dosages: seq[int8]
  ]
  VarCoef = tuple[af: float, afbin: int, coef: float]
  Model = tuple[
    id: string,
    offset: float,
    coefs: Table[string, VarCoef]   # Keyed by vid
  ]
  Models = Table[string, Model]     # Keyed by id


proc scanDosages(f: MemFile): tuple[samples: seq[string], nvariants:int] = 
  # Count lines and columns in the file.  Check for sample uniqueness.
  var
    header_read = false
    nvariants = 0
    samples = @[""]
  for slice in f.memSlices:
    if header_read == false:
      header_read = true
      samples = ($slice).strip(leading=false).split(sep="\t")[3..^1]
      var sample_set = initSet[string](sets.rightSize(samples.len))
      for sample in result.samples:
        assert sample_set.containsOrIncl(sample) == false
    else:
      nvariants += 1
  result = (samples, nvariants)


proc loadDosages(path: string, n_afbins=50): Dosages = 
  stderr.write("Loading dosages from " & path & ", " & $n_afbins & " AF bins...\n")

  # Scan the file to get sample IDs and variant count.
  # Verify sample ID uniqueness.
  let
    f = memfiles.open(path, mode=fmRead)
    (samples, nvariants) = scanDosages(f)
    nsamples = samples.len
  result.samples = samples
  stderr.write("  " & $nvariants & " variants x " & $nsamples & " samples found.  Allocating...")

  # Preallocate data
  result.dosages = newSeq[int8](nsamples*nvariants)
  result.vids = newSeq[string](nvariants)
  result.afs = newSeq[float](nvariants)
  result.afbins = newSeq[int](nvariants)
  result.vid2idx = initTable[string, int](initialSize=tables.rightSize(nvariants))
  result.afbin2idx = newSeqWith(n_afbins, newSeq[int]())
  stderr.write("  Reading...\n")

  # Read the data and verify vid uniqueness.
  var 
    i = 0
    vid_set = initSet[string](sets.rightSize(nvariants))
    header_read = false
    buffer: TaintedString = ""
  
  for line in f.lines(buffer):
    if header_read == false:
      header_read = true
      continue

    if i %% 10000 == 0:
      stderr.eraseLine()
      stderr.write("    " & $i & " / " & $nvariants & " variants")
    
    let
      fields = line.strip(leading=false).split(sep="\t")
      vid = fields[0]
    assert vid_set.containsOrIncl(vid) == false
    result.vid2idx[vid] = i
    result.vids[i] = vid
    result.afs[i] = fields[2].parseFloat
    result.afbins[i] = min(n_afbins - 1, floor(result.afs[i] * n_afbins.float).int)
    result.afbin2idx[result.afbins[i]].add(i)
    let dosages_offset = i*nsamples
    for j in 3..<fields.len:
      result.dosages[dosages_offset + (j-3)] = fields[j].parseInt.int8
    i += 1

  # Check AF bin occupancy
  var min_occupancy = result.afbin2idx[0].len
  for i in 1..<n_afbins:
    min_occupancy = min(min_occupancy, result.afbin2idx[i].len)

  stderr.write("\n  Loaded " & $nvariants & " variants x " & $nsamples & " samples, smallest AF bin size: " & $min_occupancy & "\n")


proc loadModels(path: string, n_afbins: int): Models = 
  stderr.write("Loading models from " & path & ", " & $n_afbins & " AF bins...")
  result = initTable[string, Model]()
  let f = system.open(path, mode=fmRead)
  discard f.readLine()
  for line in f.lines:
    let
      fields = line.strip(leading=false).split(sep="\t")
      model_id = fields[0]
      vid = fields[1]
      coef = fields[4].parseFloat

    if vid == "OFFSET":
      result[model_id].offset = coef
      continue

    let
      af = fields[3].parseFloat
      afbin = min(n_afbins - 1, floor(af * n_afbins.float).int)
    if not result.hasKey(model_id):
      result[model_id] = (id:model_id, offset:0.0, coefs:initTable[string, VarCoef]())

    result[model_id].coefs[vid] = (af:af, afbin:afbin, coef:coef)

  stderr.write("  Loaded " & $result.len & " models\n")


proc calcValues(model: Model, dosages: Dosages): seq[float] = 
  let n_samples = dosages.samples.len
  result = newSeq[float](n_samples)

  for i in 0..<n_samples:
    result[i] = model.offset

  for vid, coef in model.coefs.pairs:
    if dosages.vid2idx.hasKey(vid):
      let dosages_offset = dosages.vid2idx[vid]*n_samples
      for i in 0..<n_samples:
        if dosages.dosages[dosages_offset + i] == -1:
          # TODO: Can choose coef.af or dosages.afs here -- population bias vs genotype sampling bias.  GTS bias basically informative missingness.
          # coef.af solution (population bias):
          # result[i] += 2.0*coef.af*coef.coef
          # dosages.afs solution (informative missingness issues):
          result[i] += 2.0*dosages.afs[dosages.vid2idx[vid]]*coef.coef
        else:
          result[i] += float(dosages.dosages[dosages_offset + i])*coef.coef
    else:
      for i in 0..<n_samples:
        # Note use of coef.af here -- no recourse to population freq
        # because none is available.  Effect will be translation across
        # all individuals, so relative differences will be preserved.
        result[i] += 2.0*coef.af*coef.coef


proc generateNullModel(model: Model, dosages: Dosages, seed: int): Model = 
  randomize(seed)

  result.id = model.id
  result.offset = model.offset
  result.coefs = initTable[string, VarCoef](initialSize=tables.rightSize(model.coefs.len))

  for vid, coef in model.coefs.pairs:
    if dosages.vid2idx.hasKey(vid):
      # Select a new variant from dosages with matching allele frequency
      # Alternative code for first line: afbin = coef.afbin

      # Note: sampling with replacement.  Shouldn't matter in almost all cases.
      # TODO: Ideally should not select variants in LD.  Difficult to implement though.
      # One rough approach could be to enforce a minimum distance.
      # This sampling with replacement *will* be an issue for WGP.

      # TODO: Consider sampling only within xsomes.
      let
        afbin = dosages.afbins[dosages.vid2idx[vid]]
        new_vid_idx = dosages.afbin2idx[afbin][random(dosages.afbin2idx[afbin].len)]
        new_vid = dosages.vids[new_vid_idx]
      # TODO: Add random sign flipping -- good idea?  Consider the implicit bias that's
      # introduced by the "not reference" coding of the dosages.  Though OTOH this
      # should also be accompanied by AF flipping... in which case will the effect be
      # the same?  Think on this.
      result.coefs[new_vid] = (af:dosages.afs[new_vid_idx], afbin:afbin, coef:coef.coef)
    else:
      # This locus was absent from dosages, so its null equivalent should
      # be missing too.  Easily done by leaving it alone.
      discard


proc emitHeader(destination: File) = 
  destination.write("model\tsample\titer\tseed\tnafbins\tvalue\n")


proc emitValues(model: Model, dosages: Dosages, iter: int, seed: int, n_afbins: int, values: seq[float], destination: File) = 
  for i in 0..<values.len:
    destination.write(model.id & "\t" & dosages.samples[i] & "\t" & $iter & "\t" & $seed & "\t" & $n_afbins & "\t" & $values[i] & "\n")


proc calculationLoop(dosage_path: string, model_path: string, output_file: File, iters: int, n_afbins: int, seed: int) = 
  let models = loadModels(model_path, n_afbins)
  let dosages = loadDosages(dosage_path, n_afbins)

  stderr.write("Preparing random seed vector...\n")
  randomize(seed)
  var subseeds: seq[int] = @[]
  for i in 0..<iters:
    subseeds.add(random(int.high))

  stderr.write("Writing output...\n")
  emitHeader(output_file)

  var j = 0
  for model_id, model in models.pairs:
    j += 1
    stderr.eraseLine()
    stderr.write("Model " & $j & " / " & $models.len & ": " & model_id)
    let native_values = calcValues(model, dosages)
    emitValues(model, dosages, 0, seed, n_afbins, native_values, output_file)
    for i in 1..iters:
      let null_model = generateNullModel(model, dosages, subseeds[i-1])
      let null_values = calcValues(null_model, dosages)
      emitValues(model, dosages, i, seed, n_afbins, null_values, output_file)

  stderr.eraseLine()
  stderr.write("Done.\n")


proc printUsage(message: string = "") =
  if message != "":
    stderr.write(message & "\n\n")
  stderr.write("""
popstar: Calculate polygenic models and permuted nulls.

Usage: popstar [options] --dosages|d=DOSAGES --models|m=MODELS

Required parameters:
  --dosages|d=DOSAGES  Path to the input allele dosages file
  --models|m=MODELS    Path to the input model coefficient file

Options:
  --out|o=OUT     Path to the output file [default: stdout]
  --iter|i=ITER   Number of resampled null iterations to calculate [default: 1000]
  --bins|b=BINS   Number of allele frequency bins for null allele matching [default: 50]
  --seed|r=SEED   PRNG seed [default: 314159265]

v""" & VERSION & "\n" & AUTHOR & "\n\n")


proc main() =
  var
    dosage_path: string = nil
    model_path: string = nil
    output_path: string = nil
    seed: int = 314159265
    iters: int = 1000
    bins: int = 50

  for kind, key, val in getopt():
    case kind
    of cmdArgument:
      printUsage("ERROR: does not accept positional arguments.")
      return
    of cmdLongOption, cmdShortOption:
      case key
      of "dosages", "d": dosage_path = val
      of "models", "m": model_path = val
      of "out", "o": output_path = val
      of "seed", "r": seed = val.parseInt
      of "iter", "i": iters = val.parseInt
      of "bins", "b": bins = val.parseInt
      else:
        printUsage("ERROR: unrecognised option " & key)
        return
    else: assert false


  if dosage_path == nil:
    printUsage("ERROR: Dosage file path is required.")
    return

  if model_path == nil:
    printUsage("ERROR: Model file path is required.")
    return

  let output_file = if output_path == nil: stdout else: system.open(output_path, fmWrite)

  calculationLoop(dosage_path, model_path, output_file, iters, bins, seed)


main()
