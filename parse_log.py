import fileinput
import argparse
import glob
import os
from threading import Thread
import csv

def get_options():
    """
        collect options from the user
        :return: list of command line arguments
        :rtype: list
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="path to mixcr log files")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="output path for storing alignment/assemble stats")
    parser.add_argument("-pa", "--prefix_align", type=str, required=False,
                        help="prefix of alignment log file",
                        default="LOG_ALIGN")
    parser.add_argument("-a","--alignment", type=str, required=False,
                        help="alignment stats output file name",
                        default="align_stats.csv")
    parser.add_argument("-pm", "--prefix_assemble", type=str, required=False,
                        help="prefix of assemble log file",
                        default="LOG_ASSEMBLE")
    parser.add_argument("-m","--assemble", type=str, required=False,
                        help="assemble stats output file name",
                        default="assemble_stats.csv")
    return parser.parse_args()

def write_dict_to_csv(csv_file, csv_columns, dict_data):
    """
        write alignment/assemble stats to a csv file
        :param csv_file: csv file name
        :param csv_columns: header of csv file
        :param dict_data: alignment/assemble stats
        :type csv_file: str
        :type csv_columns: list
        :type dict_data: dict
    """
    try:
        with open(csv_file, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in dict_data:
                writer.writerow(data)
    except IOError as e:
            print("I/O error: {0}".format(e.message))

def process_alignment_stats(log_path, prefix, file_name):
    """
        parse alignment log file to collect alignment stats
        :param log_path: log file directory
        :param prefix: prefix of log file
        :param file_name: output file name
        :type log_path: str
        :type prefix: str
        :type file_name: str
    """
    align_files = glob.glob(os.path.join(log_path, "%s*.txt"%prefix))
    stats_rows = []
    with fileinput.input(files=align_files) as f:
        row = {}
        for line in f:
            try:
                if line.startswith("Analysis Date") or \
                    line.startswith("Overlapped") or \
                    line.startswith("Analysis time") or \
                    line.startswith("Command line arguments") or \
                    line.startswith("Paired-end alignment conflicts eliminated") or \
                    line.startswith("J gene chimeras") or \
                    line.startswith("V gene chimeras") or \
                    line.startswith("Chimeras"): continue
                if line.startswith("==="):
                    chain_list = ['TRA chains', 'TRB chains', 'TRD chains', 'TRG chains','TRA,TRD chains',
                                    'IGH chains','IGK chains', 'IGL chains']
                    for chain in chain_list:
                        if chain not in row.keys():
                            row[chain] = '0'
                    if "chains" not in row.keys():
                        row["chains"] = '0'
                    print ("align:%s" %sorted(row.keys()))
                    stats_rows.append(row)
                    row = {}
                    continue

                if not "Version" in line:
                    parts = line.strip().split(":")
                    key, val = parts[0], parts[1]
                else:
                    parts = line.strip().split(";")
                    key, val = "Version", parts[0].split(":")[1] + "," + parts[-1].split("=")[1]
                if key == 'Input file(s)':
                    files = [os.path.basename(f) for f in val.split(",")]
                    val = ",".join(files)
                elif key == 'Output file':
                    val = os.path.basename(val)
                elif "%)" in val:
                    val = val.split()[0]
                row[key] = val
            except Exception as e:
                print ("error in parsing log file: {0}".format(e.message))
                continue

    write_dict_to_csv(file_name, stats_rows[0].keys(), stats_rows)

def process_assemble_stats(log_path, prefix, file_name):
    """
        parse assemble log file to collect assemble stats
        :param log_path: log file directory
        :param prefix: prefix of log file
        :param file_name: output file name
        :type log_path: str
        :type prefix: str
        :type file_name: str
    """
    assemble_files = glob.glob(os.path.join(log_path, "%s*.txt"%prefix))
    stats_rows = []
    with fileinput.input(files=assemble_files) as f:
        row = {}
        for line in f:
            try:
                if line.startswith("Analysis Date") or \
                    line.startswith("Analysis time") or \
                    line.startswith("Command line arguments"): continue
                if line.startswith("==="):
                    chain_list = ["IGH chains", "IGK chains", "IGL chains",
                                  "TRA chains", "TRB chains", "TRD chains", "TRG chains"]
                    for chain in chain_list:
                        if chain not in row.keys():
                            row[chain] = "0"
                    print ("assemble: %s" %sorted(row.keys()))
                    stats_rows.append(row)
                    row = {}
                    continue
                if "Version" not in line:
                    parts = line.strip().split(":")
                    key, val = parts[0], parts[1]
                else:
                    parts = line.strip().split(";")
                    key, val = "Version", parts[0].split(":")[1] + "," + parts[-1].split("=")[1]
                if key == 'Input file(s)':
                    files = [os.path.basename(f) for f in val.split(",")]
                    val = ",".join(files)
                elif key == 'Output file':
                    val = os.path.basename(val)
                elif key == 'Version':
                    val = val.split(";")[0]
                elif "%)" in val:
                    val = val.split()[0]
                row[key] = val
            except Exception as e:
                print ("error in parsing log file: {0}".format(e.message))
                continue

    write_dict_to_csv(file_name, stats_rows[0].keys(), stats_rows)


def main():
    args = get_options()
    print (args)
    #start alignment thread
    print ("processing alignment stats....")
    align_thread = Thread(target=process_alignment_stats,
                          args=[args.input, args.prefix_align, os.path.join(args.output, args.alignment)])
    align_thread.start()
    #start assemble thread
    print("processing assemble stats....")
    assemble_thread = Thread(target=process_assemble_stats,
                             args=[args.input, args.prefix_assemble, os.path.join(args.output, args.assemble)])
    assemble_thread.start()

if __name__ == "__main__":
    main()
