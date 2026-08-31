#!/usr/bin/env python3
import argparse
import pysam

from tqdm import tqdm


def main(args):
  print("adding trinucleotide context")
  donor_mutation_file = add_tri_nucleotide_context(args.genomePath, args.inputPath, args.outputPath)


def add_tri_nucleotide_context(genome_file, inputPath, outputPath):
  ref = pysam.FastaFile(genome_file)
  correct = 0

  with open(inputPath, "r") as inputFile, open(outputPath, "w") as outputFile:
    i = 0
    for line in tqdm(inputFile.readlines()):

      if i > 0:
        line_list = line.split("\t")
        chrom = line_list[0]
        start = line_list[1]
        mutated_from = line_list[3]
        mutated_to = line_list[4]

        try:
          tri_ref = ref.fetch(str(chrom), int(start) - 1, int(start) + 2)
        except KeyError:
          tri_ref = ref.fetch(f"chr{str(chrom)}", int(start) - 1, int(start) + 2)
#        except ValueError:
#          tri_ref = ""

        output_line = line.split("\t")
        output_line[-1] = output_line[-1].replace("\n", "")

        if len(tri_ref) > 0:
          output_line.extend([ tri_ref[0] + mutated_from + tri_ref[2], tri_ref[0] + mutated_to + tri_ref[2], tri_ref[1] + "\n"])
        else:
          output_line.extend(["", "", "\n"]) # to fix different number of columns than expected downstream

        line = "\t".join(output_line)
        outputFile.write(line)

      else:
        output_line = line.split("\t")
        output_line[-1] = output_line[-1].replace("\n", "")
        output_line.extend(['tri_context_ref', 'tri_context_alt', 'reference_genome\n'])
        line = "\t".join(output_line)
        outputFile.write(line) # write header

      i += 1

  ref.close()


if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-g", "--genomePath", type=str, required=True, help="Help goes here") 
  parser.add_argument("-i", "--inputPath", type=str, required=True, help="Help goes here")
  parser.add_argument("-o", "--outputPath", type=str, required=True, help="Help goes here")
  args = parser.parse_args()
  main(args)
