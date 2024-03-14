import re
import argparse
import matplotlib.pyplot as plt

# Define a function to extract segments from a line
def extract_segments(line):
    # Define a regular expression pattern to match the numbers
    pattern = r'total anchors: (\d+), total segs: (\d+),.*long: (\d+) : (\d+)'

    # Use re.search to find the match in the text
    match = re.search(pattern, line)

    if match:
        total_anchors = int(match.group(1))
        total_segs = int(match.group(2))
        long_segs = int(match.group(3))
        long_anchors = int(match.group(4))

        # print("Total anchors:", total_anchors)
        # print("Total segs:", total_segs)
        # print("Long segs:", long_segs)
        # print("Long anchors:", long_anchors)
    else:
        print("Pattern not found in the text.")
        
    return total_anchors, total_segs, long_segs, long_anchors

def extract_dataset(line):
    # Define a regular expression pattern to match the values
    pattern = r'(\d+k)to(\d+k)'

    # Use re.search to find the match in the text
    match = re.search(pattern, line)

    if match:
        # Extract the matched values
        start_value = match.group(1)
        end_value = match.group(2)
        
        # print("Start value:", start_value)
        # print("End value:", end_value)
    else:
        print("Pattern not found in the text.")
    return start_value, end_value

# Create an argument parser to get the output file name from the command line
parser = argparse.ArgumentParser(description='Compute and plot a histogram of segments from an output file.')
parser.add_argument('output_file', help='Path to the output file containing segment data')
# parser.add_argument('profile_file', help='Path to the profile file containing kernel runtime')
args = parser.parse_args()

# Read the output file specified in the command line argument
tasks = {}
with open(args.output_file, 'r') as file:
    batches = []
    for line in file:
        if line.startswith("[M::main] CMD: "):
            # finish last batch 
            start_value, end_value = extract_dataset(line)
            print(f"Finish dataset {start_value} to {end_value}")
            if len(batches) > 0:
                batch_anchors = 0
                batch_long_anchors = 0
                batch_segs = 0
                batch_long_segs = 0
                for batch in batches:
                    batch_anchors += batch[0]
                    batch_segs += batch[1]
                    batch_long_segs += batch[2]
                    batch_long_anchors += batch[3]
                    print(f"Long anchor pct: {batch[3] / batch[0] * 100:.2f}%")
                    print(f"Long seg pct: {batch[2] / batch[1] * 100:.5f}%")
                print(f"Batch long anchor pct: {batch_long_anchors / batch_anchors * 100:.2f}%")
                print(f"Batch long seg pct: {batch_long_segs / batch_segs * 100:.5f}%")
                tasks[(start_value, end_value)] = batches
        if line.startswith("[M::main::"):
            print("Start new dataset...")
            batches = []
        if line.startswith("[DEBUG] total anchors"):
            total_anchors, total_segs, long_segs, long_anchors = extract_segments(line)
            batches.append((total_anchors, total_segs, long_segs, long_anchors))

# visualize the data
# plot the bar chart of the long anchor percentage and long segment percentage with respect to the dataset name
ax, fig = plt.subplots()
x = []
y = []
ticks = []

values = []
for key, value in tasks.items():
    values.append((int(key[0][:-1]), f"{key[0]}to{key[1]}", value[0][3] / value[0][0] * 100))
# sort ticks by the first element of the tuple
ticks.sort(key=lambda tup: tup[0])
for value in values:
    ticks.append(value[1])
    x.append(value[0])
    y.append(value[2])

plt.bar(x, y)
plt.xticks(x, ticks, rotation=90)
plt.xlabel("Dataset")
plt.ylabel("Long anchor percentage [%]")
plt.title("Long anchor percentage vs. dataset")
plt.tight_layout()
# save the figure to a file with the same name as the input file
output_filename = args.output_file.split('.')[0]
print(f"Saving figure to {output_filename}_long_anchor_percentage.png")
plt.savefig(f'{output_filename}_long_anchor_percentage.png')



# # Save the figure with an appropriate name based on the input file name
# output_filename = args.output_file #.split('.')[0]  # Remove the file extension
# plt.savefig(f'{output_filename}_segment_histogram.png')

# Display the histogram
# plt.show()
