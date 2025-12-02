# Run this on an HPC

## Open a terminal and allocate a node 
salloc

#Install dependencies if needed

pip install biopython primer3-py pandas matplotlib seaborn scipy

conda install -c bioconda mafft

#Clone seroBA repo

git clone https://github.com/GlobalPneumoSeq/seroba.git

cd seroba

## within the seroBA folder, save the files pangraph_primer_design.py and visualize_cps_coverage.py

module load Python/3.12.3-GCCcore-13.3.0

# pangenome approach
python pangraph_primer_design.py \
  --database /gpfs/gibbs/project/weinberger_daniel/dmw63/artic_serotype/seroba/database \
  --output ./pangraph_output

exit


# Now pare down to a small number of amplicons; looking for 5-10 amplicons per serotype that will provide redundancy and allow for discrimination

salloc --time=04:30:00

python multiplex_primer_optimizer_claude2.py \
  --database /gpfs/gibbs/project/weinberger_daniel/dmw63/artic_serotype/seroba/database \
  --primers ./pangraph_output/primers/pangraph_primers.csv \
  --output ./final_optimized_v2 \
  --max-dg -9.0 \
  --max-hairpin-tm 50.0 \
  --min-per-serotype 5
  
## Visualize coverage  
  salloc --time=04:30:00

  python visualize_cps_coverage.py \
  --database /gpfs/gibbs/project/weinberger_daniel/dmw63/artic_serotype/seroba/database \
  --primers final_optimized_v2/final_primers/optimized_multiplex_primers.csv \
  --output ./final_coverage_check \
  --plot-serotypes 19F,19A,6A,6B,37,3,7F,14,23F

## Remove overlappying amplicons

## Default (70% overlap threshold)
python deduplicate_primers.py \
  --database /gpfs/gibbs/project/weinberger_daniel/dmw63/artic_serotype/seroba/database \
  --primers final_optimized_v2/final_primers/optimized_multiplex_primers.csv \
  --output ./deduplicated_primers

## Visualize de-deuplicated amplicons

python visualize_cps_coverage.py \
  --database /gpfs/gibbs/project/weinberger_daniel/dmw63/artic_serotype/seroba/database \
  --primers deduplicated_primers/deduplicated_primers.csv \
  --output ./final_coverage_check_deduplicated \
  --plot-serotypes 19F,19A,6A,6B,37,3,7F,14,23F
  
  
  
  
  
  
  
  
  