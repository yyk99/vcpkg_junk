#!/usr/bin/env ruby

# Simple Hello World
puts "Hello, World!"

# With variables
name = "Unity Developer"
puts "Hello, #{name}!"

# Multiple lines
puts <<~MESSAGE
  Welcome to Ruby!
  This is a multi-line string.
MESSAGE

# With methods
def greet(name)
  "Hello, #{name}!"
end

puts greet("Ruby")
