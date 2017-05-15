# Alineamiento de Secuencias en Biología Molecular Computacional

El presente trabajo muestra la aplicación de algoritmos de alineación de secuencias conocidos como needleman-wunsch (global), smith-waterman (local) y semi-global con sus variantes (kband o  afín de costo por gap). Esto es muy útil en labores de biología molecular, pero requieren de alta capacidad de memoria y tiempo para ser ejecutados.

## Algoritmos de alineamiento

### Needleman-Wunsch

Para el alineamiento global de dos secuencias se empleó el algoritmo Needleman-Wunsch el cual se basa en la programación dinámica. Su funcionamiento se basa en la construcción de una matriz  F donde cada uno de sus elementos F [i, j] será el valor para el alineamiento óptimo de X [0…i] y Y [0…j]; donde X y Y corresponde a dos cadenas de texto. Para iniciar se debe de inicializar el valor de la matriz F [0, 0] con 0 y cada celda se llena con el valor óptimo de alineamiento según sea:
```
X[i] con Y[j]
X[i] con un gap
Y[j] con un gap
```
Cuando se tiene la matriz llena se procede a leer el valor de la esquina inferior derecha F [n, m] y se realiza el trace-back que consiste en agregar un carácter u otro según la ruta óptima. Dónde:
 
- Diagonal: agregar un carácter a cada hilera (a la izquierda).
- Arriba: agregar un gap a la hilera 1 y un carácter a la izquierda de la hilera 2.
- Izquierda: agregar un gap a la hilera 2 y un carácter a la izquierda de la hilera 1.

### Smith-Waterman

El alineamiento local se resuelve mediante el algoritmo Smith-Waterman, el cual consiste en una variante del algoritmo de alineamiento global. Los cambios realizados con respecto al alineamiento global son:

1. La primera columna y la primera fila inician con ceros.
2. El traceback dará inicio en la casilla con mayor valor en toda la matriz y se extenderá hasta llegar a una casilla con valor cero.

### Semi-global

El algoritmo del alineamiento semiglobal se basa en el algoritmo del alineamiento global.Este busca penalizar los gaps internos del alineamiento, no así los que se encuentre al inicio o final de alguna de las secuencias. Para ello, se debe realizar:
1. La colocación de ceros la primera fila y columna. 
2. Iniciar el traceback en el mayor valor de la última columna o fila hasta llegar a una casilla de valor cero.

### KBand

El algoritmo del alineamiento semiglobal se basa en el algoritmo del alineamiento global.Este busca penalizar los gaps internos del alineamiento, no así los que se encuentre al inicio o final de alguna de las secuencias. Para ello, se debe realizar:
1. La colocación de ceros la primera fila y columna. 
2. Iniciar el traceback en el mayor valor de la última columna o fila hasta llegar a una casilla de valor cero.