find_package (GTK2 COMPONENTS gtk)
if (GTK2_FOUND)
  message ("*** found GTK2 (enabling demos)")
  include_directories (${GTK2_INCLUDE_DIRS})
  add_executable (demo
    demo.cpp
    interaction_handle.cpp
    interactive_compound_objectives.cpp
    planar_distance.cpp
    planar_orientation_objective.cpp
    planar_robot.cpp
    )
  target_link_libraries (demo kinematic_objectives ${GTK2_LIBRARIES})
else (GTK2_FOUND)
  message (WARNING "*** GTK2 development files are required for the demos.")
endif (GTK2_FOUND)
