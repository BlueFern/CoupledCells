cmake_minimum_required(VERSION 2.8.0 FATAL_ERROR)

project(CoupledCells)

set(PROJECT_NAME "Coupled Arterial Cells")
set(PROJECT_COPYRIGHT "UC HPC, University of Canterbury, NZ")

# TODO: Update version numbers.
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 9)
set(PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
set(PROJECT_DESCRIPTION "Large-scale simulation of chemical species transport in coupled arterial cells.")

# Add local CMake directory to CMake module path to load project-specific modules. 
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

# Take care of the non-standard location of malloc includes on Mac.
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	# Location of malloc.h on OS X. 
	include_directories(/usr/include/malloc/)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Unable to run without MPI.
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})

# Unable to run without HDF5.
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIR})

# ODE implementation choice.
set(ODE_SOLVER "RK_suite" CACHE STRING
	"ODE solver implementation to be used in the simulations.")

# Available options for ODE implementation.
set(ODE_SOLVER_CHOICES "RK_suite;SUNDIALS_arkode;BOOST_odeint" CACHE INTERNAL
	"List of values for choosing a solver to be compiled against.")

# Enable ODE choices in the GUI drop-down menue.
# The appropriate choice will have to be typed in when using CLI.  
set_property(CACHE ODE_SOLVER PROPERTY STRINGS ${ODE_SOLVER_CHOICES})

# In this case we want to compile RK Suite in the contrib directory.
if(${ODE_SOLVER} STREQUAL "RK_suite")
	message(STATUS "Selected ODE solver: Local RK Suite")
	
	# Pass preprocessor macro definition to the compiler.
	add_definitions(-DRK_SUITE)
	
	# Add RK Suite code from contrib.
	include_directories(${PROJECT_SOURCE_DIR}/contrib)
	add_subdirectory(contrib)
	
	# Remember the library for the ODE solver implementation to link in to the executable.
	set(ODE_SOVLER_LIBS RKSuite)

# In this case we want to find the system SUNDIALS libraries, including ARKODE.
elseif(${ODE_SOLVER} STREQUAL "SUNDIALS_arkode")
	message(STATUS "Selected ODE solver: system SUNDIALS arkode")
	
	# Passe preprocessor macro definition to the compiler.
	add_definitions(-DARK_ODE)
	
	# Find SUNDIALS
	find_package(SUNDIALS REQUIRED)
	# Add the location of the SUNDIALS headers.
	include_directories(${SUNDIALS_INCLUDES})
	
	# Remember the libraries for the ODE solver implementation to link in to the executable.
	set(ODE_SOLVER_LIBS ${SUNDIALS_LIBRARIES})
	
# In this case we want to find the system BOOST libraries, including ODEINT.
elseif(${ODE_SOLVER} STREQUAL "BOOST_odeint")

#  _____ ___________ _____     _____ _                          _     _   _     _            _                               __   _   _            _____ ___  ___      _          _  __       _        _                              _        __                                 _                                    _      _                  _ _   _       _   _           ______                 _   _                 _      _       _                   _       _           _                   _             _     
# |_   _|  _  |  _  \  _  |_  /  ___| |                        | |   | | | |   (_)          | |                             / _| | | | |          /  __ \|  \/  |     | |        (_)/ _|     | |      | |                            (_)      / _|                               | |                                  | |    | |                (_) | | |     | | | |          | ___ \               | | ( )               | |    (_)     | |                 | |     | |         | |                 | |           | |    
#   | | | | | | | | | | | (_) \ `--.| |_ _____      ____ _ _ __| |_  | |_| |__  _ ___    ___| | __ _ _   _ ___  ___    ___ | |_  | |_| |__   ___  | /  \/| .  . | __ _| | _____   _| |_   ___| |_ __ _| |_ ___ _ __ ___   ___ _ __    _ ___  | |_ ___  _ __   _   _  ___  _   _  | |_ ___     ___ ___  _ __ ___  _ __ | | ___| |_ ___  __      ___| |_| |__   | |_| |__   ___  | |_/ / ___   ___  ___| |_|/ ___    ___   __| | ___ _ _ __ | |_ ______ _ __ ___| | __ _| |_ ___  __| |   ___ ___  _ __ | |_ ___ _ __ | |_   
#   | | | | | | | | | | | |    `--. \ __/ _ \ \ /\ / / _` | '__| __| | __| '_ \| / __|  / __| |/ _` | | | / __|/ _ \  / _ \|  _| | __| '_ \ / _ \ | |    | |\/| |/ _` | |/ / _ \ | |  _| / __| __/ _` | __/ _ \ '_ ` _ \ / _ \ '_ \  | / __| |  _/ _ \| '__| | | | |/ _ \| | | | | __/ _ \   / __/ _ \| '_ ` _ \| '_ \| |/ _ \ __/ _ \ \ \ /\ / / | __| '_ \  | __| '_ \ / _ \ | ___ \/ _ \ / _ \/ __| __| / __|  / _ \ / _` |/ _ \ | '_ \| __|______| '__/ _ \ |/ _` | __/ _ \/ _` |  / __/ _ \| '_ \| __/ _ \ '_ \| __|  
#   | | \ \_/ / |/ /\ \_/ /_  /\__/ / ||  __/\ V  V / (_| | |  | |_  | |_| | | | \__ \ | (__| | (_| | |_| \__ \  __/ | (_) | |   | |_| | | |  __/ | \__/\| |  | | (_| |   <  __/ | | |   \__ \ || (_| | ||  __/ | | | | |  __/ | | | | \__ \ | || (_) | |    | |_| | (_) | |_| | | || (_) | | (_| (_) | | | | | | |_) | |  __/ ||  __/  \ V  V /| | |_| | | | | |_| | | |  __/ | |_/ / (_) | (_) \__ \ |_  \__ \ | (_) | (_| |  __/ | | | | |_       | | |  __/ | (_| | ||  __/ (_| | | (_| (_) | | | | ||  __/ | | | |_ _ 
#   \_/  \___/|___/  \___/(_) \____/ \__\___| \_/\_/ \__,_|_|   \__|  \__|_| |_|_|___/  \___|_|\__,_|\__,_|___/\___|  \___/|_|    \__|_| |_|\___|  \____/\_|  |_/\__,_|_|\_\___| |_|_|   |___/\__\__,_|\__\___|_| |_| |_|\___|_| |_| |_|___/ |_| \___/|_|     \__, |\___/ \__,_|  \__\___/   \___\___/|_| |_| |_| .__/|_|\___|\__\___|   \_/\_/ |_|\__|_| |_|  \__|_| |_|\___| \____/ \___/ \___/|___/\__| |___/  \___/ \__,_|\___|_|_| |_|\__|      |_|  \___|_|\__,_|\__\___|\__,_|  \___\___/|_| |_|\__\___|_| |_|\__(_)
#                                                                                                                                                                                                                                                              __/ |                                            | |                                                                                                                                                                                                        
#                                                                                                                                                                                                                                                             |___/                                             |_|                                                                                                                                                                                                        
# ______           _                                  _       
# | ___ \         | |                                | |      
# | |_/ / ___  ___| |_   _ __ ___  __ _  __ _ _ __ __| |___   
# | ___ \/ _ \/ __| __| | '__/ _ \/ _` |/ _` | '__/ _` / __|  
# | |_/ /  __/\__ \ |_  | | |  __/ (_| | (_| | | | (_| \__ \_ 
# \____/ \___||___/\__| |_|  \___|\__, |\__,_|_|  \__,_|___( )
#                                  __/ |                   |/ 
#                                 |___/                       
#
# mmmmmmmdddmmmmdddhhhhhhhhhhhhhhhyyysssssssssssssssssoo++///////::::::::::-------------------...--------......-::/+ooosss
# dddddddddddddddddhhhhhhhhhhhhhhyyyysssssssssssssssssoo++/////////:::::::::------------------.....-----.......-://+ooosss
# ddddddddddddddddhhhyyyyyyhhhhhyyyssssssssoooosssssooo+++////////////::::::------------------....------.......-://+ooosss
# dddddddddddddddhhhyyyyyyyyhhhhyysssoosooooooooooooooo++//////////////:::::------------------...-------.......-://++oooos
# dddddddddddddddhhhyyyyyyyyyyyyysssoooooooooooooooooo+++///////////////::::-------------------..-------.......-://++oooos
# dddddhhhhhhhhhhhhyyyyyyyyyyyyyyssooooooooooooooooooo++++//+++o+++//////:::------------------......---........-:/++ooosss
# dddddhhyyyyyyyyyyssssyyyyyyyyysssoooooooooooooooosyhdddddhhddhyysooo+oo+::------------------.................-:/++ooosss
# ddddhyso+//+++++++++oooosssssssooo++++++++o+ossyyhhhhyyyyhhyso++//+osyhhyo+::--------------..................-//++ooosss
# hhhyso/:......-------:::://+++++////////+oyhhhhhhhysyhdmmNmmdddys++//oyddddho:-------------..................://++ooosss
# yso+//-..````````````````...-----------/shdhssyhddmmNMMMMMMMMMMMMNNmho/+syhdho:-----------...................://++ooosss
# hy+/--.......````````````````````````./yhy+/odNNMMMMMMMMMMMMMMMMMMMMMMmy+:+yhdy/:---------...................://++ooosss
# mhs/::://////:-..```````````````````:shyo/odNNMMMMMMMMMMMMMMMMMMMMMMMMMMNhooyhdh+::---::::--................-::/++ooosss
# mhyoooooossso+/:-.``.....`````````./yss+/hNNNNNMMMMMNNNNNNNNNNNNMMMMMMMMMMNmsshdho/:::/++//:-............--::///+++ooooo
# hysoooossyyysso/:...------..``````-ys++smNNNNNNNNNNNNNNNmmmmmmNNNNNNNMMMMMMMNdyhhh+::/+oo++:-....--::--.-:/+++++++++++++
# ysoooosssyyysso/:..--::/::-..`..-+hyosdNNNmmmmmmddmmmmmmmmmmNNNNNNNNNNNNMMMMMNmhhhh/::/++/:-.``.-:://::::/+ossso++//////
# hhyyyyyyyyyyyso/-..-::////:-.-:+shhsodNNmmmddhhdddmmmNNNNNNNNNNNNNNNNNNNNNMMMNNmdhhs---::--..```.-:////+oossyysoo+/:::-:
# dddddddhhhhyyso+:--://++++/:::+yhhsoymmmmmdddddmmNNNNNNNNNNNNNNNNNNNMMMMNNNNNNNmdhyy:::--:::--.``.---:+syhhhhhyso+:-----
# dddddddhyyyyyyysso+++ooooo+::/oyhyooydddddddmmmmmmmmddmmdddhhhddddhddmmNMMNNNNNmdysys+++/+oo++/-..```./shdmmddhys+:.....
# +osyyyyyyyyhhhhhhyysssssso+//+shho+oyhhhhhdmmmmhhyysssyyyyyyso+/:::/+oosydmmNNNmdysshhysoooooo+/-..```-+yddmddhyo/:.````
# ..:/osyyyyhhhdddddhhhyyssoooyhhdh/:+yhhhhdmmyo+:--..-::+shhs/:::::+ossssssyh++++hhsoyddhso+////:--..```-/oyyyyso/:-.````
# `.-/osyhhhddddddmmmdddhhyyyhddmmd::/shyyyyo+/////:------ohho-:------:/+//ooysmms/osoyddhs/-....-...``````.::/:-.````````
# `.-/+syhhhhhhddddmmmmddhhhhddmmNm:-//+++yo/:::///:-.....hNNh/:-..-.`.:oso/+ssddhy/:/ohhyo/.```````````  ````````        
# ...-:+oosssyyyyyyhhhhhhhyyhddmmmd--:`+o+o/::---.```..../mMMNy+o::/+osshhdhyyhmmNmmo+hdso/-.``````````      `````        
# ```....----:::::///++++++osyhhddh:...os+s+/:/+++//+/::+smNMMNdhyo++/++oyhdmmmNNNNmsydsy+:.```````````   ``````````      
# ```````````````````......-/+oyyhy/.`-osyssssso++////+sosdNMMNmmmmmmdmmmmmmNNNNNNNmyyd+yo:```   ```````  ```.-:::.``     
# .```.```````````.```````...-/ossyhs/-+syhdddddhhyhhhhyshdNNNMNNdhhddmmmmmNNNNNmmmmhdddsy/.```........`````.:+oo+/-``    
# .`````````````......``......-:+os//+:/syyhhdddmmmddyo+yhssysyddmh//+oyyhhdddmmmmmmmmsyhh:..-:////:::--....:+ossso/-.``` 
# ``````````````````````........-:/o/o+osysyhhhhhhyyo/::o:.-::+:-/ys+::/+oosyyyhddmmdmsydy---://++++++osyyssoossyyso+:-.``
# ``````````````````````..-...```..++o+yosssyssoooo+/:/+/-...-++oyddhs+::/+oossshdmmdmhdms-..-://///shmmNNmmdysssyyyys+:.`
# ```````````````````........``    -o+oyoooo+///++/::/ooo/-...-/oyhdhho+:/+osssyhdmmmmNdd+-.````.../hmNMMMMNmdyssyyhhys+-`
# `````````````````..-:////:-.```   /o+os+o+/:///::-:/+o+/:---:+oyyhhyy/:+syyhhhhddmmNNNyo/-``````.odNMMMMMMNmdhyyyhhhyso/
# ..`````......```..:+syhhhyso/-.```.shysooo/:/+//:---:/:---------------/+shdhhhhdddyshhyso/.`````.odNNMMMMMMNmdhysssyhhdh
# .............---:+syhdddddhys+:...`-shhs+o+/:/+++::-....-::///syhhhhys+osyhhhhhdddo:+ssso/-.````./ydmNNNMMMMNNdhso+shdmm
# ```````.....-:/osyhddmmmmmdhy+:...../hy+/++/::/++//+oooo++////+oyhddddyoosyyyyhhhh-.:+++/:-......-:+shdmNMMMMNmdy+:+shdm
# ```````````.-/oshdddmmmmmmdhy+:....`.....//////////osoo+/::-.-:+syhdddhysssssyyyho``..---....---...-:/ohmNNNNNmdho//osyh
# ````````````-/oyyhhdddmmmmdhs+-``````````.:::://://ossoo+///+oyhhhhddhyssssossyos.`````````.:///:://+++oydmNNNmmdyo++oss
# ``````````..-/osyyyyhhhhhhyso:.```````````----:::::+oossssssyyyhddhhhyyso+/++++os`````````./+ooooossyyyyydmmmmmdhyo++++/
# ````````..-:/+ossssoossssso/:.````````````-/..------/+++o+ooooosssoooo++:----:+yy```....`.:/osssyyhhddddddmmmddhso/::-..
# `````.-:/+osyyhhhhyyso++//:-.`````````````-o:.`...``.------::-::::-----....--+ydh``.--..`.-/ossyhdddddddmmmmdhs+:---....
# ```.:/+oosyhdddmmmddhyo/:::--............``o+:-.``````````````````````....-:+yhdh` ``.--.`.:/oyhddmmmddddddhy+:-..-----:
# ---/+oooossyhddmmmmmdhs+/osyyo/-------.``  +o+:-.````      ```````````.--:/+shhdh `   `....-:+oyhhhhhyssyss+:-..-/osssss
# -:/+++////+osyhhhdhhhyo+sdmmmmh+-----.`    :o+/:--.`````````````````..::+++oyhhdy  `    ```-:::/+++oo+//::-.``.:oyhhdddd
# -:/++/:::/+ooosyyyysso++ymNMMNmy:--.`      `oo+/::-....````````.....-:/+oooshhddy  `     ```-:----:++o+:.````-+yhdddddmm
# .-/++///+osyssyyyyhyyysoydNMMNho:.`         -oo+//::---......-----:::+oosssyyhdds` `       ``--``.-::/o+:.`.:ohdmmdddddd
# `-:/++++oosssssyyhhhhyyoohmdy+:.`            .+++///::--------::::://+ossyyyyhhy+`         ` `-.`  `...-:::/sydddddhyyss
# `.-:/+++++/++osyyyyysso+/oo:.``               `-/////::-------:::///+oosyyyyyss+-`            ...    `..`.:+yhhhhyyso+//
#  `.-:::::--::+osyyyso+//:-..`                   `-://::::::------:/+ooossssss+-`              `.``    `..`.-+syyso+/:---
#   ``.......--:/+osoo+/:--..`                      `.-::::------.--://++ooo+:.                  ..`      `. `--+o+::-....
#    ```..---::---:///::-.```                          `..---........--:///.`                    ````      `. `..---..````
#     ````-:/++/-......-.```                                 ``````````...`                       ``        ``` `--```````
# ````````.-//:-.`````..` `                                                                       ```        `.` `-`  ````
# ..`````.`.`````   ```.` `                                                                        ``         `` ``.   `  
# -.````````        ``.`                                                                           ``          `` `..     
# ..`````           `.`                                                                            ````         `` `.`    
# ```              ```                                                                              ``           `  ..    
#                  ``                                                                               ``            `  `.   
#                 ```                                                                                `            `` ``   
#                   `                                                                                              ````   
#                  `                                                                                                ` ``  
#                ` `                                                                                                   `  
#
# ______  _____    _____ _             ___   _____ _____ _____ _____                    _             _ _           _                      _             _                               _       _          _   _                       _                  _ _     _       _     _       _     _ _       _     _           _   
# | ___ \/  ___|  |_   _| |           / _ \ /  ___/  __ \_   _|_   _|                  | |           (_) |         | |                    | |           | |                             | |     | |        | | | |                     | |                (_) |   (_)     | |   (_)     | |   | (_)     | |   | |         | |  
# | |_/ /\ `--.     | | | |__   ___  / /_\ \\ `--.| /  \/ | |   | |    _ __   ___  _ __| |_ _ __ __ _ _| |_    __ _| |__   _____   _____  | | ___   ___ | | _____   _ __ ___  _   _  ___| |__   | |__   ___| |_| |_ ___ _ __  __      _| |__   ___ _ __    _| |_   _ ___  | |__  _  __ _| |__ | |_  __ _| |__ | |_ ___  __| |  
# |  __/  `--. \    | | | '_ \ / _ \ |  _  | `--. \ |     | |   | |   | '_ \ / _ \| '__| __| '__/ _` | | __|  / _` | '_ \ / _ \ \ / / _ \ | |/ _ \ / _ \| |/ / __| | '_ ` _ \| | | |/ __| '_ \  | '_ \ / _ \ __| __/ _ \ '__| \ \ /\ / / '_ \ / _ \ '_ \  | | __| | / __| | '_ \| |/ _` | '_ \| | |/ _` | '_ \| __/ _ \/ _` |  
# | |    /\__/ /    | | | | | |  __/ | | | |/\__/ / \__/\_| |_ _| |_  | |_) | (_) | |  | |_| | | (_| | | |_  | (_| | |_) | (_) \ V /  __/ | | (_) | (_) |   <\__ \ | | | | | | |_| | (__| | | | | |_) |  __/ |_| ||  __/ |     \ V  V /| | | |  __/ | | | | | |_  | \__ \ | | | | | (_| | | | | | | (_| | | | | ||  __/ (_| |_ 
# \_|    \____(_)   \_/ |_| |_|\___| \_| |_/\____/ \____/\___/ \___/  | .__/ \___/|_|   \__|_|  \__,_|_|\__|  \__,_|_.__/ \___/ \_/ \___| |_|\___/ \___/|_|\_\___/ |_| |_| |_|\__,_|\___|_| |_| |_.__/ \___|\__|\__\___|_|      \_/\_/ |_| |_|\___|_| |_| |_|\__| |_|___/ |_| |_|_|\__, |_| |_|_|_|\__, |_| |_|\__\___|\__,_(_)
#                                                                     | |                                                                                                                                                                                                           __/ |           __/ |                      
#                                                                     |_|                                                                                                                                                                                                          |___/           |___/                       
#

endif(${ODE_SOLVER} STREQUAL "RK_suite")

file(GLOB HDRs src/*.h)
file(GLOB SRCs src/*.cpp)

add_executable(coupledCellsModel ${SRCs} ${HDRs})

# Can use either set_target_properties or set.
# set_target_properties(coupledCellsModel PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")

# Can use either set_target_properties or set.
# set_target_properties(coupledCellsModel PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")

target_link_libraries(coupledCellsModel 
	${MPI_LIBRARIES}
	${HDF5_LIBRARIES}
	${ODE_SOVLER_LIBS})

find_package(Doxygen)

find_Package(Graphviz)
if(DOXYGEN_DOT_EXECUTABLE)
	set(HAVE_DOT "YES")
	set(CALL_GRAPH "YES")
	set(CALLER_GRAPH "YES")
else(DOXYGEN_DOT_EXECUTABLE)
	set(HAVE_DOT "NO")
	set(CALL_GRAPH "NO")
	set(CALLER_GRAPH "NO")
endif(DOXYGEN_DOT_EXECUTABLE)

set(DOXY_INPUT_DIR ${PROJECT_SOURCE_DIR}/src)
set(DOXY_OUTPUT_DIR ${PROJECT_BINARY_DIR}/doc/doxygen)
file(MAKE_DIRECTORY ${DOXY_OUTPUT_DIR})

set(DOXY_API_LINK_TAGS ${DOXY_OUTPUT_DIR}/DoxyLink.tag)

configure_file(
	${PROJECT_SOURCE_DIR}/doc/doxygen/Doxyfile.in
	${PROJECT_BINARY_DIR}/Doxyfile
	@ONLY)

if(DOXYGEN_FOUND)
	add_custom_target(doxyDocs
	COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
	COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

find_package(Sphinx)

if(NOT DEFINED SPHINX_THEME)
	set(SPHINX_THEME default)
endif()
if(NOT DEFINED SPHINX_THEME_DIR)
	set(SPHINX_THEME_DIR)
endif()

set(SPHINX_BUILD_DIR ${PROJECT_BINARY_DIR}/_sphinxbuild)
# Muffle the warning about missing directory.
# Might have to use _static directory in source to copy files here if necessary.
file(MAKE_DIRECTORY ${SPHINX_BUILD_DIR}/_static)
set(SPHINX_OUTPUT_DIR ${PROJECT_BINARY_DIR}/doc/sphinx)

configure_file(
    ${PROJECT_SOURCE_DIR}/doc/sphinx/conf.py.in
    ${SPHINX_BUILD_DIR}/conf.py
    @ONLY)

if(SPHINX_FOUND)
	add_custom_target(sphinxDocs
    COMMAND ${SPHINX_EXECUTABLE}
        -q -b html
        -c ${SPHINX_BUILD_DIR}
        -d ${SPHINX_BUILD_DIR}
        ${PROJECT_SOURCE_DIR}/doc/sphinx
        ${SPHINX_OUTPUT_DIR}
    COMMENT "Building HTML documentation with Sphinx" VERBATIM)
endif(SPHINX_FOUND)

