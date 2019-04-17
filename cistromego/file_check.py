import sys
import re
import operator
from corelib import *

### Define a check bed function here
def checkBed(file_path):
    """
    :param path of bed file:
    :return: a dictionary indicate whether a bed file is a valid and five-column bed
    """
    Info("Start checking bed file format.")
    is_5col_bed = lambda a_line: re.search(r"^chr\S+\s\d+\s\d+\s\S+\s\d+[.]*[\d]*",a_line)
    is_3col_bed = lambda a_line: re.search(r"^chr\S+\s\d+\s\d+",a_line)
    is_3tab_delimite = lambda a_line: len(a_line.split('\t')) >= 3 and len(a_line.split('\t')) < 5
    is_5tab_delimite = lambda a_line: len(a_line.split('\t')) >= 5

    def checkHeadTail(aline):
        check_res = {"fivecol": False, "pass": False}
        if is_5tab_delimite(aline) and is_5col_bed(aline):
            check_res["fivecol"] = True
            check_res["pass"] = True
        elif is_3tab_delimite(aline) and is_3col_bed(aline):
            check_res["pass"] = True
        else:
            pass
        return check_res

    with open(file_path, "rU") as pfile:
        lines = pfile.readlines()
        for line in lines:
            if line.startswith("#"):
                continue
            else:
                first_line = line
                last_line = lines[len(lines)/2]
                break
    res = {"fivecol": False, "pass": False}
    res_head = checkHeadTail(first_line)
    res_tail = checkHeadTail(last_line)
    res["fivecol"] = res_head["fivecol"] and res_tail["fivecol"]
    res["pass"] = res_head["pass"] and res_tail["pass"]
    if res["fivecol"]:
        Info("Five column bed file input.")
    else:
        Info("Your input bed file is not valid five column. Program will use all your peaks.")
    if res["pass"]:
        Info("Valid bed file format.")
    else:
        Info("ERROR: Invalid bed file format.")
        sys.exit(1)

    return res


def readPeaks(file_path, peakn):
    Info("Start reading peaks bed file.")
    check_bed_res = checkBed(file_path)
    peaks = []  # chr summit peak None
    peaki = 0
    if check_bed_res["fivecol"]:
        with open(file_path, "rU") as pfile:
            for line in pfile:
                if line.startswith("#"):
                    continue
                else:
                    peaks.append((line.rstrip().split("\t")[0],
                              (int(line.rstrip().split("\t")[1]) + int(line.rstrip().split("\t")[2])) / 2.0, 0,
                              "peak_%s" % peaki, float(line.rstrip().split("\t")[4])))
                    peaki += 1
        peaks.sort(key=lambda x: x[-1], reverse=True)
        peaks = [i[:-1] for i in peaks]
        if peakn == "all":
            peaks_top = peaks[:]
        elif float(peakn) > 1.0:
            peaks_top = peaks[:int(peakn)]
        else:
            peaks_top = peaks[:float(peakn)*peaki]

    else:
        with open(file_path, "rU") as pfile:
            for line in pfile:
                if line.startswith("#"):
                    continue
                else:
                    try:
                        peaks.append((line.rstrip().split("\t")[0],
                                  (int(line.rstrip().split("\t")[1]) + int(line.rstrip().split("\t")[2])) / 2.0, 0,
                                  "peak_%s" % peaki))
                        peaki += 1
                    except:
                        continue
        peaks_top = peaks[:]
    Info("Finish reading peaks bed file. %s peaks in total. " % peaki)
    return peaks, peaks_top


def checkExpr(file_path, expr_info):
    Info("Start checking expression file format.")
    is_refseq = lambda x: re.search(r'[A-Z][A-Z]_\d+\.*\d*(_at)?', x)
    is_genesymbol = lambda x: re.search(r'\w+', x)
    is_value = lambda x: re.search(r'\"*(\-*)(\d+)\.*(\d*)\"*', x)

    gene_column = int(expr_info.split(",")[0]) - 1
    logfc_column = int(expr_info.split(",")[1]) - 1
    fdr_column = int(expr_info.split(",")[2]) - 1

    for sep in ["\t", ",", " "]:
        res_check = {"pass": True, "sep": sep}
        with open(file_path, "rU") as exprfile:
            lines = exprfile.readlines()
            first_line = lines[1].strip().split(sep)
            last_line = lines[len(lines)/2].strip().split(sep)
        check_first_line = (is_refseq(first_line[gene_column]) or is_genesymbol(first_line[gene_column])) and is_value(first_line[logfc_column]) and is_value(first_line[fdr_column])
        check_last_line = (is_refseq(last_line[gene_column]) or is_genesymbol(last_line[gene_column])) and is_value(last_line[logfc_column]) and is_value(last_line[fdr_column])

        if check_first_line or check_last_line:
            Info("Valid expression file format. Separator is " + sep + ".")
            return res_check

    Info("ERROR: Invalid expression file format.")
    sys.exit(1)


### NOTE: add a function to translate gene ID, return the dict(geneID, symbol)
def transGeneID(expr_list, assembly):
    specie = assembly[:2]
    geneids_all = readJson('%s_geneIDs.json' % specie)
    genes_list = [i[0] for i in expr_list]
    genes_list_correct = set(map(lambda x: x.replace("_at", "").split(".")[0], genes_list))
    genes_list_length = float(len(genes_list))
    geneids_intersect = {}
    for i in geneids_all.keys():
        geneids_intersect[i] = len(set(geneids_all[i]) & genes_list_correct) / genes_list_length
    if max(geneids_intersect.values()) <= 0.6:
        Info("WARNING: Less than 50 percent of genes in the differential expression file can be mapped to gene symbol.")
        # wns += "Less than 50 percent of genes in the differential expression file can be mapped to gene symbol. "
        ### NOTE: get some warnning for users
    genesid_match = max(geneids_intersect.iteritems(), key=operator.itemgetter(1))[0]
    expr_list.sort(key=lambda x:x[-1], reverse=True)

    expr_dict_symbol = {}
    if genesid_match == "symbol":
        for i in expr_list:
            expr_dict_symbol[i[0]] = i
        return expr_dict_symbol
    else:
        genesid_match_dict = readJson('data/%s_%s_symbol.json' % (specie, genesid_match))
        for i in expr_list:
            try:
                tmp_symbol = genesid_match_dict[i[0].replace("_at", "").split(".")[0]]
                expr_dict_symbol[tmp_symbol] = [tmp_symbol] + i[1:]
            except:
                continue
        return expr_dict_symbol


### NOTE: Here to choose to seperate genes and translate the gene ID here
def readExpr(expr_path, expr_info, assembly):
    Info("Start reading expression file. ")
    check_expr_res = checkExpr(expr_path, expr_info)
    gene_column = int(expr_info.split(",")[0]) - 1
    logfc_column = int(expr_info.split(",")[1]) - 1
    fdr_column = int(expr_info.split(",")[2]) - 1
    expr_list = []
    with open(expr_path, "rU") as exprfile:
        for line in exprfile:
            try:
                tmp = line.strip().split(check_expr_res["sep"])
                expr_list.append([tmp[gene_column], float(tmp[logfc_column].replace("NA", "0")), float(tmp[fdr_column].replace("NA", "1"))])
            except:
                continue
    # print expr_list
    expr_dict_symbol = transGeneID(expr_list, assembly)
    # NOTE: add the situation equals zero
    # if len(expr_dict_symbol.values()) < 10000:
    #     Info("WARNING: Less than 10,000 genes in the differential expression gene list. ")
        # wns += "Less than 10,000 genes in the differential expression gene list. "
    # for i in expr_list_symbol:
        # if dego == "allgenes":
        # expr_dict[i[0]] = i # [geneID, logFC, FDR]
        # elif dego == "upgenes": # NOTE: change code here
        #     if i[1] > 0:
        #         expr_dict[i[0]] = [i[0], i[-1]]
        #     else:
        #         expr_dict[i[0]] = [i[0], 1.0]
        # else:
        #     if i[1] < 0:
        #         expr_dict[i[0]] = [i[0], i[-1]]
        #     else:
        #         expr_dict[i[0]] = [i[0], 1.0]
    upgenes = [i for i in expr_dict_symbol.values() if i[1]>=0]
    downgenes = [i for i in expr_dict_symbol.values() if i[1]<=0]
    Info("Finish reading expression file. %s genes in total. %s upregulated genes. %s downregulated genes. "
         % (len(expr_dict_symbol.values()), len(upgenes), len(downgenes)))
    return expr_dict_symbol, upgenes, downgenes
