using Combinatorics
using CSV
using DataFrames
using MySQL

cnx = DBInterface.connect(MySQL.Connection, "localhost", "root", "root", db = "phalp_plos", unix_socket = "/Applications/MAMP/tmp/mysql/mysql.sock")

### retrieve Phyla
data = DBInterface.execute(cnx, "SELECT DISTINCT(lineage_phylum) FROM hosts WHERE lineage_phylum != '';") |> DataFrame
phyla = data[!, 1]

### retrieve Classes
data = DBInterface.execute(cnx, "SELECT DISTINCT(lineage_class) FROM hosts WHERE lineage_class != '';") |> DataFrame
classes = data[!, 1]

### retrieve Orders
data = DBInterface.execute(cnx, "SELECT DISTINCT(lineage_order) FROM hosts WHERE lineage_order != '';") |> DataFrame
orders = data[!, 1]

### retrieve Families
data = DBInterface.execute(cnx, "SELECT DISTINCT(lineage_family) FROM hosts WHERE lineage_family != '';") |> DataFrame
families = data[!, 1]

### retrieve Genera
data = DBInterface.execute(cnx, "SELECT DISTINCT(lineage_genus) FROM hosts WHERE lineage_genus != '';") |> DataFrame
genera = data[!, 1]

### retrieve Species
data = DBInterface.execute(cnx, "SELECT DISTINCT(lineage_species) FROM hosts WHERE lineage_species != '';") |> DataFrame
species = data[!, 1]

### retrieve all architectures
all_architectures = CSV.read("/Users/stefftaelman/Downloads/Documents_local/PROJECTS/PAPER2019/simplified_architectures_condensed.csv", header = false)
endolysin_accs = Vector(CSV.read("/Users/stefftaelman/Downloads/Documents_local/PROJECTS/PAPER2019/75certain_endolysin_accs.csv", header = false).Column1)
endolysin_architectures = all_architectures[in.(all_architectures.Column1, Ref(endolysin_accs)), :]

function get_archs(cnx, host_clade::String)
    if host_clade in phyla
        if host_clade == "Actinobacteria"
            res = input("Do you mean the phylum or the class? (p/c)");
            if res == "p"
                host_string = string(" AND h.lineage_phylum = 'Actinobacteria'")
            elseif res == "c"
                host_string = string(" AND h.lineage_class = 'Actinobacteria'")
            else
                println("Please answer with either 'p' or 'c'.")
            end
        else
            host_string = string(" AND h.lineage_phylum = '", host_clade, "'")
        end
    elseif host_clade in classes
        host_string = string(" AND h.lineage_class = '", host_clade, "'")
    elseif host_clade in orders
        host_string = string(" AND h.lineage_order = '", host_clade, "'")
    elseif host_clade in families
        host_string = string(" AND h.lineage_family = '", host_clade, "'")
    elseif host_clade in genera
        host_string = string(" AND h.lineage_genus = '", host_clade, "'")
    elseif host_clade in species
        host_string = string(" AND h.lineage_species = '", host_clade, "'")
    elseif host_clade == "Cyanophyceae"
        host_string = string("AND h.lineage_phylum = 'Cyanobacteria'")
    elseif host_clade == "all"
        host_string = ""
    else
        println("Please enter a valid host clade to filter on.")
    end
    #=
    if protein_type == "all"
        type_string = ""
    elseif protein_type == "endolysin"
        type_string = "up.type = 'endolysin' AND "
    elseif protein_type == "VAPGH"
        type_string = "up.type = 'VAPGH' AND "
    else
        println("Please enter a valid type to filter on: all, endolysin or VAPGH")
    end
    =#
    ### get accessions
    #query = string("SELECT DISTINCT(up.UniProt_ID) FROM UniProt as up JOIN link_phage_host as l JOIN hosts as h WHERE ", type_string, "up.phages_ID = l.phages_ID AND l.hosts_ID = h.hosts_ID", host_string, ";")
    query = string("SELECT DISTINCT(up.UniProt_ID) FROM UniProt as up JOIN link_phage_host as l JOIN hosts as h WHERE up.phages_ID = l.phages_ID AND l.hosts_ID = h.hosts_ID", host_string, ";")
    data = DBInterface.execute(cnx, query) |> DataFrame
    accs = data[!, :UniProt_ID]

    ### get subsetted dataframe
    archs_df = endolysin_architectures[in(accs).(endolysin_architectures.Column1), :]

    ### get architectures
    arch = Vector{}()
    for i in 1:size(archs_df)[1]
        doms = [k for k in collect(archs_df[i, 2:end]) if !(ismissing(k))]
        if length(doms) > 1
            row = ""
            for (jdx, j) in enumerate(doms)
                row = string(row, j, ";")
                if jdx == length(doms)
                    row = [row[1:end-1]]
                end
            end
        else
            row = doms
        end
        append!(arch, row)
    end
    return arch #.|> a -> Symbol.(split(a, ";"))
end

function get_subset_properties(architectures)
    ### get subset properties
    lens = Vector{Integer}()
    domains = Set{String}()
    for i in architectures
        append!(lens, length(split(i, ";")))
        [push!(domains, j) for j in split(i, ";")]
    end
    modules = collect(domains)
    max_len = maximum(lens)
    #return println("Branch contains architectures composed from $(length(modules)) domains in a maximum of $(max_len) positions.")
    return modules, max_len
end

function symbolize(archs, mod2sym, sym2mod)
    moduledb = Vector{String}()
    for j in archs
        if length(j) == 1
            append!(moduledb, string(j))
        else
            line = ""
            for k in split(j, ";")
                line = line*string(mod2sym[k])
            end
            push!(moduledb, string(line))
        end
    end
    return moduledb
end

function iterate_regex_2_wild(modules, f::Function;
                modsperslot=2,
                verbose=false,
                kwargs...)
    best_score = -Inf
    best_reg = r""
    wildcards = ["", ".{0,1}", ".{0,2}", ".{0,3}"]
    for n1 in 1:modsperslot
        for c1 in combinations(modules, n1)
            for s1 in wildcards
                for n2 in 1:modsperslot
                    for c2 in combinations(modules, n2)
                        reg = Regex("[" * prod(c1) * "]")
                        f_reg = f(reg; kwargs...)
                        if f_reg > best_score
                            best_score = f_reg
                            best_reg = reg
                            verbose && println("$best_reg (score=$best_score)")
                        end
                        reg = Regex("[" * prod(c1) * "]" * s1 * "[" * prod(c2) * "]")
                        f_reg = f(reg; kwargs...)
                        if f_reg > best_score
                            best_score = f_reg
                            best_reg = reg
                            verbose && println("$best_reg (score=$best_score)")
                        end
                    end
                end
            end
        end
    end
    return best_reg, best_score
end

function iterate_regex_3_tame(modules, f::Function;
                modsperslot=2,
                verbose=false,
                kwargs...)
    best_score = -Inf
    best_reg = r""
    for n1 in 1:modsperslot
        for c1 in combinations(modules, n1)
            for n2 in 1:modsperslot
                for c2 in combinations(modules, n2)
                    for m1 in ["", "{0,1}"]
                        for n3 in 1:modsperslot
                            for c3 in combinations(modules, n2)
                                reg = Regex("[" * prod(c1) * "]")
                                f_reg = f(reg; kwargs...)
                                if f_reg > best_score
                                    best_score = f_reg
                                    best_reg = reg
                                    verbose && println("$best_reg (score=$best_score)")
                                end
                                reg = Regex("[" * prod(c1) * "]" * "[" * prod(c2) * "]" * m1)
                                f_reg = f(reg; kwargs...)
                                if f_reg > best_score
                                    best_score = f_reg
                                    best_reg = reg
                                    verbose && println("$best_reg (score=$best_score)")
                                end
                                reg = Regex("[" * prod(c1) * "]" * "[" * prod(c2) * "]" * m1 * "[" * prod(c3) * "]")
                                f_reg = f(reg; kwargs...)
                                if f_reg > best_score
                                    best_score = f_reg
                                    best_reg = reg
                                    verbose && println("$best_reg (score=$best_score)")
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return best_reg, best_score
end

function iterate_working(modules, f::Function;
                modsperslot=3,
                verbose=false,
                kwargs...)
    best_score = -Inf
    best_reg = r""
    for n1 in 1:modsperslot
        for c1 in combinations(modules, n1)
            for n2 in 1:modsperslot
                for c2 in combinations(modules, n2)
                    for m1 in ["", "{0,1}"]
                        for n3 in 1:modsperslot
                            for c3 in combinations(modules, n2)
                                reg = Regex("[" * prod(c1) * "]")
                                f_reg = f(reg; kwargs...)
                                if f_reg > best_score
                                    best_score = f_reg
                                    best_reg = reg
                                    verbose && println("$best_reg (score=$best_score)")
                                end
                                reg = Regex("[" * prod(c1) * "]" * "[" * prod(c2) * "]" * m1)
                                f_reg = f(reg; kwargs...)
                                if f_reg > best_score
                                    best_score = f_reg
                                    best_reg = reg
                                    verbose && println("$best_reg (score=$best_score)")
                                end
                                reg = Regex("[" * prod(c1) * "]" * "[" * prod(c2) * "]" * m1 * "[" * prod(c3) * "]" * m1)
                                f_reg = f(reg; kwargs...)
                                if f_reg > best_score
                                    best_score = f_reg
                                    best_reg = reg
                                    verbose && println("$best_reg (score=$best_score)")
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return best_reg, best_score
end

function scoring(regex; db=db, modules=modules, ent_sym=ent_sym, sym2mod=sym2mod)
    true_archs = unique(db) #relevant items
    true_pos = 0
    for i in true_archs
        if occursin(regex, i) == true
            true_pos += length(match(regex, i).match)/length(i)
        end
    end

    ### recall
    recall = true_pos/length(true_archs)

    ### precision
    matches = 0
    for i in ent_sym
        if occursin(regex, i) == true
            matches += length(match(regex, i).match)/length(i)
        end
    end
    precision = true_pos/matches
    F_score = 2*precision*recall/(precision+recall)
    return F_score
end

function desymbolize(regex, sym2mod)
    regex_string = string(regex)[3:end-1]
    dr = ""
    or_criterion = false
    for (idx, i) in enumerate(regex_string)
        if i in keys(sym2mod)
            if regex_string[idx-1] != '{' && regex_string[idx-1] != ','
                dr = dr * sym2mod[i]
                if or_criterion == true && regex_string[idx+1] != ']'
                    dr = dr * "|"
                end
            else
                dr = dr * string(i)
            end
        else
            dr = dr * string(i)
            if i == '['
                or_criterion = true
            elseif i == ']'
                or_criterion = false
            end
        end
    end
    return dr
end

function main_design_exhaustive_iterator(cnx, host_clades, iterator)
    ### get symbols
    entire = unique(get_archs(cnx, "all"))
    all_m, _ = get_subset_properties(entire)
    symbols = String(vcat(UInt8.(65:90), UInt8.(97:122)))
    nmodules = length(all_m)
    mod2sym = Dict{String, Char}()
    sym2mod = Dict{Char, String}()
    for i in 1:nmodules
        mod2sym[all_m[i]] = Char(symbols[i])
        sym2mod[Char(symbols[i])] = all_m[i]
    end
    ent_sym = symbolize(entire, mod2sym, sym2mod)

    for i in host_clades
        ### find rule
        archs = get_archs(cnx, i)
        db = symbolize(archs, mod2sym, sym2mod)
        modules = unique(string(db...))
        regex, score = iterator(modules, scoring, modsperslot=3, verbose=true, db=db, modules=modules, ent_sym=ent_sym, sym2mod=sym2mod)
        rule  = desymbolize(regex, sym2mod)

        ### write to file
        clade_string = ">" * i * "\n"
        open("/Users/stefftaelman/Downloads/Documents_local/PROJECTS/PAPER2019/designs_3pos_working.txt", "a") do io
            write(io, clade_string)
            write(io, (rule*"\n"))
            write(io, (string(score)*"\n"))
            write(io, (string(length(archs))*"\n\n"))
        end
    end
end

#=
list = ["Bacillus cereus", "Bacillus subtilis", "Bacillus thuringiensis", "Listeria monocytogenes", "Staphylococcus aureus", "Enterococcus faecalis",
        "Lactobacillus plantarum", "Lactococcus lactis", "Streptococcus agalactiae", "Streptococcus pneumoniae", "Streptococcus pyogenes",
        "Streptococcus thermophilus"]
=#

list = ["Paenibacillus larvae", "Paenibacillus", "Paenibacillaceae", "Streptococcus dysgalactiae", "Streptococcus suis"]

main_design_exhaustive_iterator(cnx, list, iterate_working) #started 07/12 17:03

b = get_archs(cnx, "Synechococcus")
c = get_archs(cnx, "Synechococcaceae")

b == c

c = get_archs(cnx, "Vibrio")
d = get_archs(cnx, "Vibrionaceae")

c == d
