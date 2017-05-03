from cursesmenu import *
from cursesmenu.items import *
from pairwise import algorithms

def submenu_item_sequences(parent_menu):
    submenu_2 = CursesMenu("Ingreso de Secuencias", "Elija una opción")
    seq1 = "SEND"
    seq2 = "AND"
    # TODO cambiar
    seq1 = list(seq1)
    seq2 = list(seq2)
    function_item_2 = FunctionItem("Needleman", algorithms.needleman_wunsch, [seq1,seq2])
    item2 = MenuItem("Another Item")
    submenu_2.append_item(function_item_2)
    submenu_2.append_item(item2)
    submenu_item_2 = SubmenuItem("Ingresar secuencias", submenu=submenu_2)
    submenu_item_2.set_menu(parent_menu)
    return submenu_item_2

def main():
    menu = CursesMenu("Alineamientos", "Elija una opción")
    
    
    menu.append_item(submenu_item_sequences(menu))
    menu.start()
    menu.join()


if __name__ == "__main__":
    main()