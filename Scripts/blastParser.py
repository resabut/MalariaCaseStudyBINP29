#! /usr/bin/python3
"""
blastParser.py

Description: Parses a .blastp file and outputs a tab-delimited table with an entry for each hit

User-defined functions: None
Non-standard modules: None

Procedure:
    1. Loops through the file
    2. Saves each query id until it finds a hit
    3. For each hit it records the target id, e-value, identity % and score of the match
    4. Write everything to the output file

Input: .blastp file
Output: Tab-separated hit table

Usage:
    python3 blastParser.py file.blastp [-o code output file]

Version: 1.00
Date: 2023-01-24
Name: Joan Escrivà Font
"""

import argparse
import re
from pathlib import Path

# Parses arguments from command line
parser = argparse.ArgumentParser()
# input file
parser.add_argument("blastp", type=Path, metavar="blastp_input_file",
                    help=".blastp input file separated by tabs")

# output_file file
parser.add_argument("-o", type=Path, dest="output_file",
                    required=False, metavar=".txt output_file file",
                    default=Path("output_file.txt"),
                    help="output_file file path")

args = parser.parse_args()

# MAIN CODE
print("Parsing file...")
with open(args.blastp) as blastp_input, open(args.output_file, "w") as output:
    # write header to output_file
    header = "#query\ttarget\te-value\tidentity(%)\tscore"
    output.write(header)
    output.write("\n")

    # Count number of hits and queries
    hit_num = 0
    query_num = 0

    # loops through the input file
    for line in blastp_input:

        # looks for the query
        query = re.search(r"Query=\s(\w*)\s",  # RegEx: Query= (query id)
                          line)

        # If the line contains a query
        #   saved the id and goes to the next line
        if query:
            query_id = query.group(1)
            query_num += 1
            continue

        # if the lines contains a hit
        elif line.startswith(">"):
            hit_num += 1  # counts one more hit
            hit = line.strip(">\n")  # remove extra chars

            # skip lines until score-expect line
            next(blastp_input)
            next(blastp_input)

            # get score and e-value
            line = next(blastp_input)
            score = re.search(r"Score = [\d.]* bits \((\d*)\)",  # RegEx: Score = bit-score bits (Score)
                              line).group(1)
            e_value = re.search(r"Expect = ([\w\-.]*),",  # RegEx: Expect = e-value (it can have exponential values)
                                line).group(1)

            # get identity % in the next line
            line = next(blastp_input)
            identity = re.search(r"Identities = [\d/]* \((\d*)%\)",  # RegEx: fraction (identity percentage)
                                 line).group(1)

            # write to output
            result = "\t".join([query_id, hit, e_value, identity, score])
            output.write(result)
            output.write("\n")

        # query has no result
        elif line.startswith("***"):
            # print to output empty result
            no_result_line = "\t".join([query_id, "", "", "", ""])
            output.write(no_result_line)
            output.write("\n")

print("DONE")
print(f"{hit_num} hits and {query_num} queries parsed")
