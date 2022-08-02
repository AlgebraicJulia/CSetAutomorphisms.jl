@info "Building nauty..."
nautydir = joinpath(@__DIR__, "nauty27r3")
bashit(str) = run(`bash -c "$str"`)

bashit("cd $nautydir; ./configure; make -j4")
