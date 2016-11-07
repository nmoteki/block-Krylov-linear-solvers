TARGET = block_iterative_solvers_test
CXX = g++-6
CXXFLAGS = -O3
CPPFLAGS = -DEIGEN_NO_DEBUG -DNDEBUG
INCLUDES = -I/Users/moteki/eigen_3_2_10
$(TARGET) : block_iterative_solvers_test.o bl_cocg_rq.o bl_bicg_rq.o bl_bicr_rq.o bl_bicgstab_rq.o bl_bicgstab.o bl_bicggr.o qr_reduced.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -o $@ block_iterative_solvers_test.o bl_cocg_rq.o bl_bicg_rq.o bl_bicr_rq.o bl_bicgstab_rq.o bl_bicgstab.o bl_bicggr.o qr_reduced.o

block_iterative_solvers_test.o : block_iterative_solvers_test.cpp block_iterative_solvers.hpp linear_algebra_addon.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c block_iterative_solvers_test.cpp

bl_cocg_rq.o : bl_cocg_rq.cpp linear_algebra_addon.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c bl_cocg_rq.cpp

bl_bicg_rq.o : bl_bicg_rq.cpp linear_algebra_addon.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c bl_bicg_rq.cpp

bl_bicr_rq.o : bl_bicr_rq.cpp linear_algebra_addon.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c bl_bicr_rq.cpp

bl_bicgstab_rq.o : bl_bicgstab_rq.cpp linear_algebra_addon.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c bl_bicgstab_rq.cpp

bl_bicgstab.o : bl_bicgstab.cpp linear_algebra_addon.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c bl_bicgstab.cpp

bl_bicggr.o : bl_bicggr.cpp linear_algebra_addon.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c bl_bicggr.cpp


qr_reduced.o: qr_reduced.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) -c qr_reduced.cpp

clean:
	rm -f *.o $(TARGET)
