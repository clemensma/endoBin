def coverage_estimate():
    with open('alignment_rate.txt', 'r') as alignment_rate_file:                # Get alignment rate from bowtie2 output
        line = alignment_rate_file.readlines()[-1]
    alignment_rate = ''
    for i in line:
        if i == '%': break
        alignment_rate += i
    alignment_rate = float(alignment_rate)/100

    with open('host_count.txt', 'r') as host_count:                             # Get no. of bases in reference genome
        reference_base_count = int(host_count.readline())


    with open('base_count.txt', 'r') as base_count:                             # Get no. of bases in read files
        reads_base_count = int(base_count.readline()) * 2

    coverage = (alignment_rate * reads_base_count) / reference_base_count       # Calculate coverage

    with open('coverage.txt', 'w') as coverage_file:
        coverage_file.write(str(coverage))

def main():
    coverage_estimate()

    
if __name__ == "__main__":
    main()
