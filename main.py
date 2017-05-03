from pairwise import algorithms
from Bio import AlignIO, SeqIO


def menu_principal():
    """Menú Principal.

    Muestra el menú principal de opciones.

        Retorna
        ------
        int
                Opción ingresada por el usuario.
        """
    print("-------------------Menú-------------------")
    print("[1] Carga de Datos")
    print("[2] Aplicar Needleman-Wunsch")
    print("[S] Salir")
    return input("Ingrese la opción que desea: ")


def submenu_item_sequences():
    print("-------------------Elija un método-------------------")
    print("[1] Por Archivos")
    print("[2] Tecleadas interactivamente")
    return input("Ingrese la opción que desea: ")


def caso_1():
    print("-------------------Secuencia Alfa-------------------")
    option = submenu_item_sequences()
    if(option == '1'):
        file_path = input("Ingrese el path: ")
        file_type = input("Ingrese el tipo: ")
        try:
            alpha = AlignIO.read(file_path, file_type)
        except:
            print("Inválido, vuelva a intentar...")
            return
    if(option == '2'):
        print("Interactivamente")
    print("-------------------Secuencia Beta-------------------")
    option = submenu_item_sequences()
    if(option == '1'):
        file_path = input("Ingrese el path: ")
        file_type = input("Ingrese el tipo: ")
        try:
            beta = AlignIO.read(file_path, file_type)
        except:
            print("Inválido, vuelva a intentar...")
            return
    if(option == '2'):
        print("Interactivamente")


def caso_2():
    seq1 = "SEND"
    seq2 = "AND"
    # TODO cambiar
    seq1 = list(seq1)
    seq2 = list(seq2)
    algorithms.needleman_wunsch(seq1, seq2)


def caso_3():
    print("caso 3")


def caso_salir():
    print("Fin del programa")


def caso_invalido():
    print('Opción Inválida. Intente de nuevo.')


# Main
OPCIONES = {'1': caso_1, '2': caso_2, '3': caso_3,
            's': caso_salir, 'S': caso_salir}
while True:
    OPCION = menu_principal()
    f = OPCIONES.get(OPCION, caso_invalido)
    f()
    if OPCION == 's' or OPCION == 'S':
        break  # break del programa
