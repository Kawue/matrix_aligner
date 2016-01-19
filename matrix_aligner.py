from __future__ import print_function
__author__ = 'Karsten'
from sys import argv
import numpy as np
import pandas as pd


def main(string_x, string_y, type, score_copy, cost_sub, cost_indel, cost_matrix):
    if type == 'n':
        needleman_wunsch(string_x, string_y)

    dataframe_embedding()


def needleman_wunsch(string_x, string_y,score_copy, cost_sub, cost_indel):
    len_x = len(string_x)+1
    len_y = len(string_y)+1
    matrix = np.zeros(shape=(len_x, len_y), dtype=int)

    for i in range(0,len_x):
        for j in range(0,len_y):
            if i == 0 and j > 0:
                matrix[0][j] = matrix[0][j-1]+cost_indel
            elif i > 0 and j == 0:
                matrix[i][0] = matrix[i-1][0] + cost_indel
            elif i > 0 and j > 0:
                if string_x[i-1] == string_y[j-1]:
                    matrix[i][j] = max(matrix[i-1][j-1] + score_copy, matrix[i-1][j] + cost_indel, matrix[i][j-1] + cost_indel)
                elif string_x[i-1] != string_y[j-1]:
                    matrix[i][j] = max(matrix[i-1][j-1] + cost_sub, matrix[i-1][j] + cost_indel, matrix[i][j-1] + cost_indel)

    return matrix


def free_end_gap(string_x, string_y,score_copy, cost_sub, cost_indel):
    len_x = len(string_x)+1
    len_y = len(string_y)+1
    matrix = np.zeros(shape=(len_x, len_y), dtype=int)

    for i in range(0, len_x):
        for j in range(0, len_y):
            if i == 0 and j > 0:
                matrix[0][j] = 0
            if i > 0 and j == 0:
                matrix[i][0] = 0
            elif i > 0 and j > 0:
                if string_x[i-1] == string_y[j-1]:
                    matrix[i][j] = max(matrix[i-1][j-1] + score_copy, matrix[i-1][j] + cost_indel, matrix[i][j-1] + cost_indel)
                elif string_x[i-1] != string_y[j-1]:
                    matrix[i][j] = max(matrix[i-1][j-1] + cost_sub, matrix[i-1][j] + cost_indel, matrix[i][j-1] + cost_indel)

    return matrix

def smith_waterman(string_x, string_y, score_copy, cost_sub, cost_indel):
    len_x = len(string_x)+1
    len_y = len(string_y)+1
    matrix = np.zeros(shape=(len_x, len_y), dtype=int)

    for i in range(0, len_x):
        for j in range(0, len_y):
            if i == 0 and j > 0:
                matrix[0][j] = 0
            if i > 0 and j == 0:
                matrix[i][0] = 0
            elif i > 0 and j > 0:
                if string_x[i-1] == string_y[j-1]:
                    matrix[i][j] = max(0, matrix[i-1][j-1] + score_copy, matrix[i-1][j] + cost_indel, matrix[i][j-1] + cost_indel)
                elif string_x[i-1] != string_y[j-1]:
                    matrix[i][j] = max(0, matrix[i-1][j-1] + cost_sub, matrix[i-1][j] + cost_indel, matrix[i][j-1] + cost_indel)

    return matrix

###needs some work
def sellers(string_x, string_y, score_copy, cost_sub, indel, errors):
    len_x = len(string_x)+1
    len_y = len(string_y)+1
    matrix = np.zeros(shape=(len_x, len_y), dtype=int)
    for i in range(0,len(matrix)):
        for j in range(0,len(matrix[i])):
            matrix[i][j] = -1
    max_range_x = None
    last_index = []
    for i in range(len_x):
        if i*indel <= errors:
            matrix[i][0] = i * indel
    last_index.append(int(errors/indel))
    for j in range(1,len_y):
        matrix[0][j] = 0
        max_range_x = min(len_x-1, last_index[j-1] + 1)
        #+1 because max_range_x has to be included
        for i in range(1, max_range_x+1):
            if string_x[i-1] == string_y[j-1]:
                if matrix[i][j-1] != -1:
                    matrix[i][j] = min(matrix[i-1][j-1]+score_copy, matrix[i][j-1]+indel, matrix[i-1][j]+indel)
                else:
                    matrix[i][j] = min(matrix[i-1][j-1]+score_copy, matrix[i-1][j]+indel)
            elif string_x[i-1] != string_y[j-1]:
                if matrix[i][j-1] != -1:
                    matrix[i][j] = min(matrix[i-1][j-1]+cost_sub, matrix[i][j-1]+indel, matrix[i-1][j]+indel)
                else:
                    matrix[i][j] = min(matrix[i-1][j-1]+cost_sub, matrix[i-1][j]+indel)
        ### When does this happen!?!
        if matrix[max_range_x][j] < errors:
            last_index.append(min(len_x-1, max_range_x+int((errors-matrix[max_range_x][j])/indel)))
            #+1 because you have to start behind max_range_x
            for i in range(max_range_x+1, last_index[j]):
                matrix[i][j] = matrix[i-1][j]+indel
        else:
            to_append = 0
            #+1 because max_range_x has to be included
            for i in range(0, max_range_x + 1):
                if matrix[i][j] <= errors:
                    to_append = i
                    #hier stimmt was noch nicht
            last_index.append(to_append)
    return matrix


def gotoh(string_x, string_y):
    len_x = len(string_x)+1
    len_y = len(string_y)+1
    matrix_s, matrix_v ,matrix_h = np.zeros(shape=(len_x, len_y), dtype=int), \
                                   np.zeros(shape=(len_x, len_y), dtype=int), \
                                   np.zeros(shape=(len_x, len_y), dtype=int)


    return matrix_s, matrix_v, matrix_h

print(gotoh("asdasd", "asdadsa"))


def create_latex_output(string_x, string_y, matrix ,type, path):
    len_x = len(string_x)
    len_y = len(string_y)
    f = open(path, "w")
    print(r'\textbf{' + str(type).capitalize() + r':}\\', file=f)

    if type == 'matrix':
        tex_columns = r'\begin{tabular}{'
        #String_x Spalte + Epsilon Spalte + String_y
        for i in range(len_y+2):
            tex_columns = tex_columns + '|c'
        else:
            tex_columns = tex_columns +  r'|}'
        print(str(tex_columns), file=f)
        print(r'\hline',file=f)

        tex_row_string_y = r' & $\epsilon$'
        for i in range(len_y):
            tex_row_string_y = tex_row_string_y + ' & $' + string_y[i] + '$'
        else:
            tex_row_string_y = tex_row_string_y +  r'\\'
        print(tex_row_string_y, file=f)
        print(r'\hline',file=f)

        for i in range(0,len_x+1):
            matrix_row = ''
            for j in range(len_y+1):
                matrix_row = matrix_row + ' & $'+str(matrix[i][j])+'$'
            else:
                matrix_row = matrix_row + r'\\'
            if i == 0:
                print(r'$\epsilon$' + matrix_row, file=f)
                print(r'\hline',file=f)
            else:
                print('$' + string_x[i-1] + '$' + matrix_row, file=f)
                print(r'\hline',file=f)
        else:
            print(r'\end{tabular}\\',file=f)
    #untested
    elif type == 'sellers':
        tex_columns = r'\begin{tabular}{'
        #String_x Spalte + Epsilon Spalte + String_y
        for i in range(len_y+2):
            if i == 0:
               tex_columns = tex_columns + 'c|'
            else:
                tex_columns = tex_columns + 'c'
        else:
            tex_columns = tex_columns + '}'
        print(str(tex_columns), file=f)

        tex_row_string_y = r' & $\epsilon$'
        for i in range(len_y):
            tex_row_string_y = tex_row_string_y + ' & $' + string_y[i] + '$'
        else:
            tex_row_string_y = tex_row_string_y +  r'\\'
        print(tex_row_string_y, file=f)
        print(r'\hline',file=f)

        for i in range(0,len_x+1):
            matrix_row = ''
            for j in range(len_y+1):
                if matrix[i][j] != -1:
                    matrix_row = matrix_row + ' & $'+str(matrix[i][j])+'$'
                else:
                    matrix_row = matrix_row + ' & '
            else:
                matrix_row = matrix_row + r'\\'
            if i == 0:
                print(r'$\epsilon$' + matrix_row, file=f)
            else:
                print('$' + string_x[i-1] + '$' + matrix_row, file=f)
        else:
            print(r'\end{tabular}\\',file=f)


def dataframe_embedding(matrix):
    dataframe = pd.DataFrame

    return dataframe

create_latex_output("AABB", "BABAABABB", sellers("AABB", "BABAABABB", 0, 1, 1, 1), 'sellers', 'C:\\Users\\Karsten\\Desktop\\Texkram.txt')