"""
Author: Levi Malmström
Run "include("Packages.jl") to install all necessary packages.
"""

using Pkg
Pkg.add("Interpolations")
Pkg.add("Integrals")
Pkg.add("Images")
Pkg.add("ImageMagick")
Pkg.add("ImageView")
Pkg.add("FileIO")
Pkg.add("DataFrames")
Pkg.add("CSVFiles")
Pkg.add("StaticArrays")
Pkg.add("SpecialFunctions")
Pkg.add("CUDA")

"""
Note: DataFrames and CSVFiles can be removed by commenting out "using DataFrames" and "using CSVFiles" from "initialization.jl"
and then delete everything after the line that says "#decide whether to run tests".
Also, if you can't/don't want to use the GPU mode, remove CUDA from the packages in this file, as well as in "initialization.jl".
At some point I will make this automatic or semi-automatic.
"""