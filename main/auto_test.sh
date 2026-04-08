#!/usr/bin/expect
#Author: Levi Malmstrom

set timeout 360

#Sphere in cartesian coordinates (LttM)
spawn julia -t auto
expect "julia>"
send -- include("initialization.jl")\r
expect "cpu or gpu"
send "cpu\r"
expect "automated test:"
send "yes\r"
expect "e.g."
send "cart_sph\r"
expect "julia>"
send "exit()\r"

#Sphere in cartesian coordinates (5P)
expect ":main"
spawn julia -t auto
expect "julia>"
send -- include("initialization.jl")\r
expect "cpu or gpu"
send "gpu\r"
expect "automated test:"
send "yes\r"
expect "e.g."
send "cart_sph\r"
expect "julia>"
send "exit()\r"

#Sphere in spherical coordinates (LttM)
expect ":main"
spawn julia -t auto
expect "julia>"
send -- include("initialization.jl")\r
expect "cpu or gpu"
send "cpu\r"
expect "automated test:"
send "yes\r"
expect "e.g."
send "sph_sph\r"
expect "julia>"
send "exit()\r"

#Sphere in spherical coordinates (5P)
expect ":main"
spawn julia -t auto
expect "julia>"
send -- include("initialization.jl")\r
expect "cpu or gpu"
send "gpu\r"
expect "automated test:"
send "yes\r"
expect "e.g."
send "sph_sph\r"
expect "julia>"
send "exit()\r"

#Schwarzschild BH with simple disk (LttM)
expect ":main"
spawn julia -t auto
expect "julia>"
send -- include("initialization.jl")\r
expect "cpu or gpu"
send "cpu\r"
expect "automated test:"
send "yes\r"
expect "e.g."
send "bh_simple\r"
expect "julia>"
send "exit()\r"

#Schwarzschild BH with simple disk (5P)
expect ":main"
spawn julia -t auto
expect "julia>"
send -- include("initialization.jl")\r
expect "cpu or gpu"
send "gpu\r"
expect "automated test:"
send "yes\r"
expect "e.g."
send "bh_simple\r"
expect "julia>"
send "exit()\r"

#DNEG Wormhole (LttM)
expect ":main"
spawn julia -t auto
expect "julia>"
send -- include("initialization.jl")\r
expect "cpu or gpu"
send "cpu\r"
expect "automated test:"
send "yes\r"
expect "e.g."
send "wormhole\r"
expect "julia>"
send "exit()\r"

#DNEG Wormhole (5P)
expect ":main"
spawn julia -t auto
expect "julia>"
send -- include("initialization.jl")\r
expect "cpu or gpu"
send "gpu\r"
expect "automated test:"
send "yes\r"
expect "e.g."
send "wormhole\r"
expect "julia>"
send "exit()\r"

#Kerr BH with simple disk and alpha = 0.998 M (LttM)
expect ":main"
spawn julia -t auto
expect "julia>"
send -- include("initialization.jl")\r
expect "cpu or gpu"
send "cpu\r"
expect "automated test:"
send "yes\r"
expect "e.g."
send "kerr\r"
expect "julia>"
send "exit()\r"

#Kerr BH with simple disk and alpha = 0.998 M (5P)
expect ":main"
spawn julia -t auto
expect "julia>"
send -- include("initialization.jl")\r
expect "cpu or gpu"
send "gpu\r"
expect "automated test:"
send "yes\r"
expect "e.g."
send "kerr\r"
expect "julia>"
send "exit()\r"

  expect eof
