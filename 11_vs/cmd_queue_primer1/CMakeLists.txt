cmake_minimum_required(VERSION 3.15)
project(cmd_queue_primer1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the executable
add_executable(cmd_queue_primer1 main.c)

# On Linux/BSD, sys/queue.h is typically available in system headers
# No additional libraries needed
