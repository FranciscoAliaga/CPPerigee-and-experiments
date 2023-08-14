# localiza el directorio con el código de fuente
src_dir = .\..\..\src\
# NOTA: puede ser necesario cambiar src si el experimento está en una carpeta más profunda

# opciones especiales para el compilador g++.
# -I src se encarga de incluir los archivos necesarios.
# -W** exigen que se entregue la mayor cantidad de avisos posibles, para diagnosticar fallas.
# -O3 exige el más alto nivel de optimización al compilador
compiler_settings = -fdiagnostics-color=always -I $(src_dir) -Wall -Weffc++ -Wextra -Wsign-conversion -Wpedantic -O3

# compila la simulación
simulation:
	g++ sim_edge_priority.cpp $(compiler_settings) -o simulacion_edge_priority

