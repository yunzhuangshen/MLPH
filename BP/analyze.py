from utils import *
from make_figure_curve import solving_gap_curve
from make_table import make_tables


def analyze(prefix='./'):
    method1='BP_MLPH_10._1._0.1'; method2='BP_def'
    make_tables(method1, method2, prefix)
    solving_gap_curve(method1, method2, prefix)

def analyze_MLPH_variants(prefix='./'):

    method1='BP_def'; method2='BP_none'
    make_tables(method1, method2, prefix)
    method1='BP_MLPH_10._1._0.1'; method2='BP_MLPH_force_exact'
    make_tables(method1, method2, prefix)
    method1='BP_MLPH_10._1._0.1'; method2='BP_MLPH_10._1._1.'
    make_tables(method1, method2, prefix)
    method1='BP_MLPH_10._1._0.1'; method2='BP_MLPH_1._1._0.1'
    make_tables(method1, method2, prefix)
    method1='BP_MLPH_10._1._0.1'; method2='BP_MLPH_0.1_1._0.1'
    make_tables(method1, method2, prefix)

if __name__ == '__main__':
    analyze()
    analyze_MLPH_variants()