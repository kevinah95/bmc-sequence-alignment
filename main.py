#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from algorithms.needleman_wunsch import plot_nw,needleman_wunsch

# TODO see more:
# https://github.com/biopython/biopython/blob/master/Bio/pairwise2.py
score_only = False
# Penalty scores
penalty = {'MATCH': 1, 'MISMATCH': -1, 'GAP': -2}
# List of the algorithms
algorithms = ['Needleman–Wunsch (Global)',
              'Smith–Waterman (Local)',
              'Semiglobal',
              'Global con K-Band']
# Sequences
sequence_alpha = None
sequence_beta = None

# Main definition - constants
menu_actions = {}

# =======================
#     MENUS FUNCTIONS
# =======================

# Main menu


def main_menu():
    os.system('clear')

    print ("Bienvenido,\n")
    print ("Por favor ingrese alguna opción:")
    print ("a. Ayuda")
    print ("b. Tablas")
    print ("c. Listar")
    print ("d. Val")
    print ("e. Match")
    print ("f. Mismatch")
    print ("g. Gap")    
    print ("h. Secuencia ALPHA")
    print ("i. Secuencia BETA")
    print ("j. Needleman-Wunsch")
    print ("0. Salir")
    choice = input(" >>  ")
    exec_menu(choice)

    return

# Execute menu


def exec_menu(choice):
    os.system('clear')
    ch = choice.lower()
    ch = ch.split()
    len_of_ch = len(ch)
    if len_of_ch == 0:
        menu_actions['main_menu']()
    if len_of_ch == 1:
        try:
            menu_actions[ch[0]]()
        except KeyError:
            print ("Invalid selection, please try again.\n")
            menu_actions['main_menu']()
    if len_of_ch == 2 and ch[0] == 'a':
        menu_actions[ch[0]]('nw')

    return

# =========
# a.Ayuda
# =========


def opt_ayuda(algorithm=''):
    if algorithm == '':
        help(opt_ayuda)
    if algorithm == 'needleman-wunsch':
        print("Call help NW")
    print ("9. Regresar")
    print ("0. Salir")
    choice = input(" >>  ")
    exec_menu(choice)
    return

# =========
# b.Tablas
# =========


def change_score_only():
    global score_only
    print("Desea cambiar el Bit de estado (s/n):")
    choice = input(" >>  ")
    choice = choice.lower()
    if (choice == "s"):
        score_only = not score_only
        print("Bit cambiado")
        show_score_only()


def show_score_only():
    global score_only
    print("Estado del Bit (Mostrar tablas): ", end='')
    print("No" if score_only else "Sí")


def opt_tablas():
    global score_only
    show_score_only()
    change_score_only()
    back_or_exit_choice()
    return

# =========
# c.Listar
# =========


def show_algorithms():
    print ('%s' % '\n'.join(map(str, algorithms)))


def opt_algoritmos():
    show_algorithms()
    back_or_exit_choice()
    return


def show_penalty():
    print(penalty)
    return

# =========
# d.Val
# =========


def opt_val():
    show_penalty()
    back_or_exit_choice()
    return

# =========
# e.Match
# =========


def give_me_a_number():
    while True:
        try:
            return int(input("Enter a number: "))
        except:
            print("Vuelva a intentar...")
            pass


def change_one_penalty(name):
    global penalty
    print("Desea cambiar el valor (s/n):")
    choice = input(" >>  ")
    choice = choice.lower()
    if (choice == "s"):
        value = give_me_a_number()
        penalty[name] = value
        print("Valor cambiado")
        show_one_penalty(name)


def show_one_penalty(name):
    print("Valor del " + name + ": ", penalty[name])


def opt_match():
    penalty_name = 'MATCH'
    show_one_penalty(penalty_name)
    change_one_penalty(penalty_name)
    back_or_exit_choice()
    return

# =========
# f.Mismatch
# =========


def opt_mismatch():
    penalty_name = 'MISMATCH'
    show_one_penalty(penalty_name)
    change_one_penalty(penalty_name)
    back_or_exit_choice()
    return

# =========
# g.Gap
# =========


def opt_gap():
    penalty_name = 'GAP'
    show_one_penalty(penalty_name)
    change_one_penalty(penalty_name)
    back_or_exit_choice()
    return

# =========
# Metodos para Ingresar secuencias
# =========

def method_file(obj):
    file_path = obj
    try:
        print("Ingrese el formato del archivo: ")
        file_format = input(" >>  ")
        first_record = next(SeqIO.parse(file_path, file_format))
    except:
        print("Vuelva a intentarlo")
        pass
    return first_record

def select_method_of_input(name):
    sequence = ""
    print("Para la secuencia ",name)
    print("Ingrese la secuencia o el path del archivo:")
    obj = input(" >>  ")
    if (os.path.exists(obj)):
        sequence = method_file(obj)
    else:
        obj = obj.upper().split()  # upper sequence and split to get first
        sequence = Seq(obj[0])
    return sequence

# =========
# h. Secuencia ALPHA
# =========
def show_alpha_sequence():
    global sequence_alpha
    print("Secuencia ALPHA")
    print(sequence_alpha if sequence_alpha else "No se ha ingresado la secuencia")
    print("")
    print("Desea cambiar el valor (s/n):")
    choice = input(" >>  ")
    choice = choice.lower()
    if (choice == "s"):
        sequence_alpha = select_method_of_input("ALPHA")

def opt_change_alpha_sequence():
    show_alpha_sequence()
    back_or_exit_choice()
    return

# =========
# i. Secuencia BETA
# =========

def show_beta_sequence():
    global sequence_beta
    print("Secuencia BETA")
    print(sequence_beta if sequence_beta else "No se ha ingresado la secuencia")
    print("")
    print("Desea cambiar el valor (s/n):")
    choice = input(" >>  ")
    choice = choice.lower()
    if (choice == "s"):
        sequence_beta = select_method_of_input("BETA")

def opt_change_beta_sequence():
    show_beta_sequence()
    back_or_exit_choice()
    return

# =========
# j. Needleman-Wunsch
# =========
#@profile
def opt_needleman_wunsch():
    if score_only :
        needleman_wunsch.needleman_wunsch(sequence_alpha,sequence_beta,penalty,score_only)
    else:
        plot_nw.plot_nw(sequence_alpha,sequence_beta,penalty)
    back_or_exit_choice()
    return


# =========
# Back or Exit
# =========

def back_or_exit_choice():
    print("---")
    print ("9. Regresar")
    print ("0. Salir")
    choice = input(" >>  ")
    exec_menu(choice)
    return


def back():
    menu_actions['main_menu']()

# Exit program


def exit():
    print("Realizado por:")
    print("-Kevin Hernández Rostrán")
    print("-Diana Vargas Hernández")
    print("Introducción a la Biología Molecular Computacional")
    print("Profesor: Esteban Arias Méndez")
    print("Periodo: I Semestre, 2017")
    print("Escuela de Computación, Instituto Tecnológico de Costa Rica")
    sys.exit()

# =======================
#    MENUS DEFINITIONS
# =======================


# Menu definition
menu_actions = {
    'main_menu': main_menu,
    'a': opt_ayuda,
    'b': opt_tablas,
    'c': opt_algoritmos,
    'd': opt_val,
    'e': opt_match,
    'f': opt_mismatch,
    'g': opt_gap,
    'h': opt_change_alpha_sequence,
    'i': opt_change_beta_sequence,
    'j': opt_needleman_wunsch,
    '9': back,
    '0': exit,
}

# =======================
#      MAIN PROGRAM
# =======================

# Main Program
if __name__ == "__main__":
    # Launch main menu
    main_menu()
