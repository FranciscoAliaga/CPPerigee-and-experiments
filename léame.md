
- Cada carpeta corresponde a un conjunto de experimentos con un interés específico.
	- la carpeta tendrá un archivo "contexto.md" con lo necesario para entender los experimentos realizados adentro.
	- cada conjunto tendrá subcarpetas con simulaciones, estas subcarpetas contendrán:
	- código:
		- makefile
		- sim.cpp
		- aquí, makefile es un archivo para el comando "make" que se encarga de compilar la simulación adecuadamente. 
		- después de compilar con make, habrá un ejecutable "simulacion".
	- carpetas con resultados:
		- gráfos
		- curvas
		- logs.


el comando make podría ejecutarse con el programa "make". En MSYS2:

pacman -S mingw-w64-x86_64-make

lo instala en:
C:\msys64\mingw64\bin\mingw32-make.exe