# localiza el directorio con el código de fuente
src_dir = .\..\..\src\
# NOTA: puede ser necesario cambiar src si el experimento está en una carpeta más profunda

# aquí capturamos los mensajes y errores de compilación
output_file = out.log

# opciones especiales para el compilador g++.
# -I src se encarga de incluir los archivos necesarios.
# -W** exigen que se entregue la mayor cantidad de avisos posibles, para diagnosticar fallas.
# -O3 exige el más alto nivel de optimización al compilador
compiler_settings = -fdiagnostics-color=always -I $(src_dir) -Wall -Weffc++ -Wextra -Wsign-conversion -Wpedantic -O3

# compila la simulación
simulation:
	g++ sim.cpp $(compiler_settings) -o simulacion > $(output_file)

# el ejecutable se llama "simulacion"