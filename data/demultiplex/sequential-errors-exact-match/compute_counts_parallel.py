import multiprocessing
from collections import defaultdict

reads = defaultdict(list)

# Function to process each line
def process_line(line):
    return line.split('\n')[-1]

# Function to process chunks of lines
def process_lines_chunk(lines):
    return [process_line(line) for line in lines]

# Function to read file in chunks
def chunk_file(fname, size=1024*1024):
    file_end = file_size(fname)
    with open(fname, 'r') as f:
        chunk_end = f.tell()
        while True:
            chunk_start = chunk_end
            f.seek(size, 1)
            f.readline()
            chunk_end = f.tell()
            yield chunk_start, chunk_end - chunk_start
            if chunk_end >= file_end:
                break

# Get file size
def file_size(fname):
    with open(fname, 'rb') as f:
        f.seek(0, 2)
        return f.tell()

# Function to read a specific chunk of lines from the file
def read_chunk(fname, start, size):
    with open(fname, 'r') as f:
        f.seek(start)
        lines = f.read(size).splitlines()
        return lines

if __name__ == '__main__':
    # Path to the file
    fname = 'large_file.txt'

    # Create a pool of processes
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())

    # Create tasks - each task is to process a chunk of lines
    tasks = [(fname, start, size) for start, size in chunk_file(fname)]

    # Map tasks to the processing functions
    results = pool.starmap(process_lines_chunk, tasks)

    # Combine results
    processed_lines = [line for sublist in results for line in sublist]

    # Output results
    for line in processed_lines:
        print(line)
