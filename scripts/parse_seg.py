import re
import argparse
import matplotlib.pyplot as plt

# Define a function to extract segments from a line
def extract_segments(line):
    match = re.search(r'long segments (\d+)', line)
    if match:
        if int(match.group(1)) > 1000000:
            print(line)
        return int(match.group(1))
    return None

# Define a function to extract runtime from a line
def extract_runtime(line):
    match = re.search(r'last launch runtime: (\d+\.\d+) ms', line)
    if match:
        return float(match.group(1))
    return None


# Initialize variables to store segment counts
segment_counts = []
runtimes = []

# Create an argument parser to get the output file name from the command line
parser = argparse.ArgumentParser(description='Compute and plot a histogram of segments from an output file.')
parser.add_argument('output_file', help='Path to the output file containing segment data')
parser.add_argument('runtime_file', help='Path to the file containing runtime data')
args = parser.parse_args()

# Read the output file specified in the command line argument
with open(args.output_file, 'r') as file:
    for line in file:
        segments = extract_segments(line)
        if segments is None:
            continue
        segment_counts.append(segments)
        

# Read the runtime file specified in the command line argument
with open(args.runtime_file, 'r') as runtime_file:
    for line in runtime_file:
        runtime = extract_runtime(line)
        if runtime is None:
            continue
        runtimes.append(runtime)

# Calculate the total number of segments
total_segments = sum(segment_counts)
total_runtime = sum(runtimes)
throughput = total_segments / total_runtime  # anchors/ms

# Create a histogram
plt.hist(segment_counts, bins=200, edgecolor='k')
plt.xlabel('Segments')
plt.ylabel('Frequency')
plt.title(f'Total Segments: {total_segments}, throughput: {throughput} anchors/ms')
plt.grid(True)

# Save the figure with an appropriate name based on the input file name
output_filename = args.output_file #.split('.')[0]  # Remove the file extension
plt.savefig(f'{output_filename}_segment_histogram.png')

# Display the histogram
# plt.show()
