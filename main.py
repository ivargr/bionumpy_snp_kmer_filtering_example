import time
from dataclasses import dataclass
import bionumpy as bnp
from bionumpy.bnpdataclass import bnpdataclass
from bionumpy.io.delimited_buffers import DelimitedBuffer
import npstructures as nps


def get_snp_kmers(fasta_filename, vcf_filename, k):
    genome = bnp.Genome.from_file(fasta_filename, filter_function=None)
    variants = genome.read_locations(vcf_filename, has_numeric_chromosomes=False)

    # get only snps from vcf
    snp_mask = (variants.get_data_field('ref_seq').shape[-1] == 1) & (variants.get_data_field('alt_seq').shape[-1] == 1)
    variants = variants[snp_mask]

    # Get windows around these snps
    windows = variants.get_windows(flank=k-1)

    # Use the windows to extra sequences (kmers) fr
    sequences = genome.read_sequence()[windows]
    sequences = bnp.as_encoded_array(sequences, bnp.DNAEncoding)
    reference_kmers = bnp.get_kmers(sequences, k)

    # Extract kmers with alternative allele on SNPs
    sequences[:, k-1] = variants.get_data_field('alt_seq').ravel()
    sequences = bnp.as_encoded_array(sequences, bnp.DNAEncoding)
    alt_kmers = bnp.get_kmers(sequences, k)

    return reference_kmers, alt_kmers


@bnpdataclass
class SequenceEntry:
    sequence: bnp.DNAEncoding


class TxtKmerBuffer(DelimitedBuffer):
    dataclass = SequenceEntry


if __name__ == "__main__":

    fasta_filename = "Mguttatus_256_v2.0.fa"
    vcf_filename = "real_grand_small-1.vcf"
    kmers_file_name = "TEST_kmerlist.txt"

    ref_kmers, alt_kmers = get_snp_kmers(fasta_filename, vcf_filename, k=31)

    # Create a HashSset of alt kmers, lets us quickly check whether any other kmers is in the set
    alt = nps.HashSet(alt_kmers.raw().ravel())

    # create a new file that we can write the filtered kmers to
    out_file = bnp.open("out.txt", "w", buffer_type=TxtKmerBuffer)

    # read all kmers from txt file, read them chunk by chunk to keep memory low
    chunks = bnp.open(kmers_file_name, buffer_type=TxtKmerBuffer).read_chunks()
    t = time.perf_counter()
    for chunk in chunks:
        # chunk.sequence contains all the sequences encoded using DNAEncoding,
        # we encode these to kmers
        chunk_kmers = bnp.get_kmers(chunk.sequence, k=31).raw().ravel()

        # Check how many of these are in our alt_kmers hashset
        is_in_alt_kmers = alt.contains(chunk_kmers)

        # is_in_alt_kmers is a boolean mask, True for every kmer that is in alt kmers, False for the rest
        # create a new chunk with those that are not in hashset
        filtered = chunk[~is_in_alt_kmers]
        print(f"Chunk contained {len(chunk)} kmers. After filtering there are {len(filtered)} left")
        out_file.write(filtered)

    out_file.close()
