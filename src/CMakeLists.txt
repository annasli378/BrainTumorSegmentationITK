cmake_minimum_required(VERSION 2.8)
PROJECT(POMwJO_Laborki) # nazwa solucji
FIND_PACKAGE(ITK REQUIRED) # tworzony projekt wymaga ITK, CMake ma go znaleźć
IF( ITK_FOUND )
INCLUDE (${ITK_USE_FILE}) # oraz dodać do ustawień projektu
ELSE( ITK_FOUND )
MESSAGE(FATAL_ERROR "Nie znaleziono ITK! Ustaw prawidłowo ITK_DIR")
ENDIF( ITK_FOUND )
ADD_EXECUTABLE(
POMwJO_Lab_2 # nazwa projektu oraz wynikowego pliku .exe
main.cxx # plik z implementacją main()
)
TARGET_LINK_LIBRARIES(
POMwJO_Lab_2 # nazwa projektu, taka jak w ADD_EXECUTABLE(...)
${ITK_LIBRARIES} # informacje dla linkera umożliwiające użycie ITK
)