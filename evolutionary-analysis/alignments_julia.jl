using DelimitedFiles;
using BioSequences;
using BioAlignments;
using CSV;
using DataFrames;

data = readdlm("/Users/stefftaelman/Downloads/Documents_local/PROJECTS/PAPER2019/PhaLP_unique_seqs.csv", ',', '\n')
seqs = data[:, 1]

### started around 10:05pm on 29/04
sim_matrix = zeros(length(seqs), length(seqs))
for i in 1:length(seqs)
    for j in i:length(seqs)
        tmp = pairalign(LocalAlignment(), seqs[i], seqs[j], AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1); score_only=true)
        sim_matrix[i, j] = score(tmp)
        sim_matrix[j, i] = score(tmp)
    end
end  #started 30/12 around 14:45

CSV.write("/Users/stefftaelman/Downloads/Documents_local/PROJECTS/PAPER2019/similarity_matrix_010120.csv",  DataFrame(sim_matrix), writeheader=false)
