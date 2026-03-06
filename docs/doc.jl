#!/usr/bin/env julia

"""
    Documentation Generation Script for CTDirect.jl

This script generates the documentation for CTDirect.jl and then removes
CTDirect from the docs/Project.toml to keep it clean.

Usage (from any directory):
    julia docs/doc.jl
    # OR
    julia --project=. docs/doc.jl
    # OR  
    julia --project=docs docs/doc.jl

The script will:
1. Activate the docs environment
2. Add CTDirect as a development dependency in docs environment
3. Generate the documentation using docs/make.jl
4. Remove CTDirect from docs/Project.toml
5. Clean up the docs environment
"""

using Pkg

println("🚀 Starting documentation generation for CTDirect.jl...")

# Step 0: Activate docs environment (works from any directory)
docs_dir = joinpath(@__DIR__)
println("📁 Activating docs environment at: $docs_dir")
Pkg.activate(docs_dir)

# Step 1: Add CTDirect as development dependency
println("📦 Adding CTDirect as development dependency...")
# Get the project root (parent of docs directory)
project_root = dirname(docs_dir)
Pkg.develop(path=project_root)

# Step 2: Generate documentation
println("📚 Building documentation...")
include(joinpath(docs_dir, "make.jl"))

# Step 3: Remove CTDirect from docs environment
println("🧹 Cleaning up docs environment...")
Pkg.rm("CTDirect")

println("✅ Documentation generated successfully!")
println("📖 Documentation available at: $(joinpath(docs_dir, "build", "index.html"))")
println("🗂️  CTDirect removed from docs/Project.toml")