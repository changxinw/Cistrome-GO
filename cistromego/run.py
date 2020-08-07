from rpscore import calcRPscore
from enrichment import annotIndices
from opt_validator import optValidate
import json

def runCistromeGO(args):
    args = optValidate(args)
    symbol_rp_dict = calcRPscore(args.bed, args.peakn, args.assembly, args.decay)
    json.dump(symbol_rp_dict, open("%s_rp.json"%args.prefix, "w"))
    if not args.expr:
        res = annotIndices(symbol_rp_dict, args.assembly, True, args.max_gene_number, args.min_gene_number, args.prefix)
    else:
        from expr_combine import calcRPT
        symbol_rpt_dict = calcRPT(args.expr, args.exprinfo, symbol_rp_dict, args.assembly, args.dego, args.logfc_cut, args.fdr_cut)
        res = annotIndices(symbol_rpt_dict, args.assembly, False, args.max_gene_number, args.min_gene_number, args.prefix)
    return res




