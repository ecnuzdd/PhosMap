from xml.etree.ElementTree import iselement
import xml.etree.ElementTree as ET
import os, sys, re

EXP_NAMES = sys.argv[1]
MASCOT_XML_DIR = sys.argv[2]
MASCOT_TXT_DIR = sys.argv[3]
COLNAMES = ['pep_seq', 'pep_query_num', 'query_rank', 'pep_score', 'pep_var_mod', 'pep_var_mod_pos', 'pep_var_mod_conf']
TAG_PREFIX = '{http://www.matrixscience.com/xmlns/schema/mascot_search_results_2}'
tag_queries = TAG_PREFIX + 'queries'
tag_query = TAG_PREFIX + 'query'
tag_queries_query = tag_queries + '/' + tag_query
tag_q_peptide = TAG_PREFIX + 'q_peptide'
tag_queries_query_q_peptide = tag_queries_query + '/' + tag_q_peptide

def get_mascot_xml_txt_file_names(mascot_xml_dir):
    mascot_xml_file_names = os.listdir(mascot_xml_dir)
    pattern = r'^E[\w]*_F[\d]_R[\d].xml$'
    target_mascot_xml_file_names = []
    output_mascot_txt_file_names = []
    for mascot_xml_file_name in mascot_xml_file_names:
        if re.match(pattern, mascot_xml_file_name):
            target_mascot_xml_file_names.append(mascot_xml_file_name)
            mascot_txt_file_name = mascot_xml_file_name.replace('.xml', '.txt')
            output_mascot_txt_file_names.append(mascot_txt_file_name)
    return target_mascot_xml_file_names, output_mascot_txt_file_names

def formatting_pep_seq(pep_seq, pep_var_mod_pos):
    #pep_seq = 'SGSISVK'
    pep_seq_new = ''
    one_flag = '0'
    # pep_var_mod_pos = '0.0010000.0'
    pos_flag = pep_var_mod_pos.split('.')[1]
    pos_flag_len = len(pos_flag)
    for i in range(pos_flag_len):
        tmp_pos_flag = pos_flag[i]
        tmp_pep_seq = pep_seq[i]
        if tmp_pos_flag != one_flag:
            tmp_pep_seq = tmp_pep_seq.lower()
        pep_seq_new = pep_seq_new + tmp_pep_seq
    return pep_seq_new

def parserMascotXMLAndWriteToFile(inputFileName, outputFileName):
    try:
        with open(outputFileName, 'w') as fw:
            fw.write('\t'.join(COLNAMES) + '\n')
            tree = ET.ElementTree(file=inputFileName)
            for elem in tree.iterfind(tag_queries_query_q_peptide):
                query_num = elem.attrib['query']
                query_rank = elem.attrib['rank']
                # queries -> query -> q_peptide
                tag_pep_seq = TAG_PREFIX + 'pep_seq'
                tag_pep_score = TAG_PREFIX + 'pep_score'
                tag_pep_var_mod = TAG_PREFIX + 'pep_var_mod' # Phospho
                tag_pep_var_mod_pos = TAG_PREFIX + 'pep_var_mod_pos'
                tag_pep_var_mod_conf = TAG_PREFIX + 'pep_var_mod_conf'

                e_tag_pep_seq = elem.find(tag_pep_seq)
                e_tag_pep_score = elem.find(tag_pep_score)
                e_tag_pep_var_mod = elem.find(tag_pep_var_mod)
                e_tag_pep_var_mod_pos = elem.find(tag_pep_var_mod_pos)
                e_tag_pep_var_mod_conf = elem.find(tag_pep_var_mod_conf)

                pep_seq = e_tag_pep_seq.text
                pep_score = e_tag_pep_score.text

                if iselement(e_tag_pep_var_mod) and e_tag_pep_var_mod.text is not None:
                    pep_var_mod = e_tag_pep_var_mod.text
                else:
                    pep_var_mod = ''

                if iselement(e_tag_pep_var_mod_pos) and e_tag_pep_var_mod_pos.text is not None:
                    pep_var_mod_pos = e_tag_pep_var_mod_pos.text
                else:
                    pep_var_mod_pos = ''

                if iselement(e_tag_pep_var_mod_conf) and e_tag_pep_var_mod_conf.text is not None:
                    pep_var_mod_conf = e_tag_pep_var_mod_conf.text
                else:
                    pep_var_mod_conf = ''

                if pep_var_mod_pos != '':
                    pep_seq = formatting_pep_seq(pep_seq, pep_var_mod_pos)

                lst = [pep_seq,  query_num, query_rank, pep_score, pep_var_mod, pep_var_mod_pos, pep_var_mod_conf]
                lst_str = '\t'.join(lst)
                fw.write(lst_str + '\n')

    except RuntimeError:
        return False
    else:
        return True


def parser_an_experiment(target_mascot_xml_file_names, output_mascot_txt_file_names, exp_name):
    target_mascot_xml_file_names_len = len(target_mascot_xml_file_names)
    for i in range(target_mascot_xml_file_names_len):
        inputFileName = target_mascot_xml_file_names[i]
        inputFilePath = os.path.join(MASCOT_XML_DIR, exp_name, inputFileName)
        outputFileName = output_mascot_txt_file_names[i]
        outputFilePath = os.path.join(MASCOT_TXT_DIR, exp_name, outputFileName)
        flag = parserMascotXMLAndWriteToFile(inputFilePath, outputFilePath)
        if flag:
            print(exp_name + ';' + inputFileName + ';' + '1')
        else:
            print(exp_name + ';' + inputFileName + ';' + '0')



def __main__():
    exp_names = EXP_NAMES.split(';')
    for exp_name in exp_names:
        mascot_xml_dir = os.path.join(MASCOT_XML_DIR, exp_name)
        target_mascot_xml_file_names, output_mascot_txt_file_names = get_mascot_xml_txt_file_names(mascot_xml_dir)
        mascot_txt_dir = os.path.join(MASCOT_TXT_DIR, exp_name)
        if not os.path.exists(mascot_txt_dir):
            os.mkdir(mascot_txt_dir)
        parser_an_experiment(target_mascot_xml_file_names, output_mascot_txt_file_names, exp_name)


if __name__ == "__main__":__main__()
