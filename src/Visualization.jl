
using Catlab.WiringDiagrams, Catlab.Graphs

function to_parts(xs::Vector{Int}; distinguish::Int=0)::String
  join([let xi=(findall(==(i), xs)), x=join(xi, " "), b=i==distinguish;
          b ? "<b>"*x*"</b>" : x end
        for i in 1:max0(xs)]," | ")
end

function view_cdict(cd::CDict; distinguish::Pair{Symbol, Int}=(:X=>0))::String
  ds, di = distinguish
  join(["$k:$(to_parts(v; distinguish=k==ds ? di : 0))" for (k,v) in pairs(cd)]," <BR/> ") * " "
end

"""Adds tree to wiring diagram, returns box id of root"""
function add_tree!(wd::WiringDiagram, t::Tree)::Int
  n = length(t.children)
  bx= Box(view_cdict(t.coloring), [Int],repeat([Int], n))
  root = add_box!(wd, bx)
  for (i,(k,v)) in enumerate(t.children)
    child = add_tree!(wd, v)
    add_wire!(wd, Wire(k, (root, i), (child, 1)))
  end
  return root
end

function view_tree(t::Tree)::WiringDiagram
  wd = WiringDiagram([],[])
  add_tree!(wd, t)
  return wd
end

struct VizState
  wd::WiringDiagram
  tree::Tree
  msg::String
  boxids::Dict{VPSI,Int}
  colors::Dict{Int,String}
end

VizState(wd::WiringDiagram,tree::Tree,msg::String,boxids::Dict{VPSI,Int}) = VizState(wd,tree,msg,boxids,Dict{Int,String}())


function to_vizstate(hs::Vector{History}; sat::Bool=false)::Vector{VizState}
  vs = VizState[]
  add_history!(vs, hs; sat)
  return vs
end

function add_history!(vs::Vector{VizState}, hs::Vector{History}; sat::Bool=false)
  if isempty(hs)
    return nothing
  end
  h = hs[1]
  v = isempty(vs) ? nothing : deepcopy(vs[end])

  if (h.action == "start")
    isnothing(v) || error()
  elseif h.action == "start_iter"
    # unpack history
    split_seq, color_seq, new_ind, splt_data = h.val
    split_tab, split_color, split_inds = splt_data
    # Load data from previous vizstate
    if isnothing(v)
      wd, vt, bi, c = WiringDiagram([],[]), Tree(), Dict{VPSI,Int}(), Dict{Int,String}()
    else
      wd, vt, bi, c = v.wd, v.tree, v.boxids, Dict([
        b=>(clr == "lightblue" ? "bisque" : clr) for (b,clr) in v.colors])
    end

    # Add new box
    bi[split_seq] = boxid = add_box!(wd, Box(view_cdict(color_seq[1]),[Int],[Int]))
    c[boxid] = "lightblue"

    if !isempty(split_seq)
      add_wire!(wd, Wire(split_seq[end][2], (bi[split_seq[1:end-1]], 1), (boxid, 1)))
    end

    # Update tree
    vt[split_seq].coloring = color_seq[1]
    vt[split_seq].saturated = color_seq[end]
    vt[split_seq].indicator = new_ind
    for i in split_inds
      vt[split_seq].children[split_tab=>i] = Tree()
    end

    push!(vs, VizState(deepcopy(wd), vt, "Add node", bi, c))
      # Add states for each color saturation step
    if sat
      for (i, color) in enumerate(color_seq[2:end])
        set_subpart!(wd.diagram, boxid, :value,  view_cdict(color))
        msg = "Color saturate $i/$(length(color_seq)-1)"
        push!(vs, VizState(deepcopy(wd), vt, msg, bi, deepcopy(c)))
      end
    end
    # Add state for splitting
    vt.saturated, vt.indicator = color_seq[end], new_ind
    set_subpart!(wd.diagram, boxid, :value,  view_cdict(
      color_seq[end]; distinguish=(split_tab => split_color)))
    col = deepcopy(c)
    if split_tab == :_nothing # we have a leaf
      msg = "Found leaf"
      col[boxid] = "green"
    else
      msg = "Split on $splt_data"
    end
    push!(vs, VizState(deepcopy(wd), vt, msg, bi, col))
  elseif h.action == "auto_prune"
    push!(vs, VizState(v.wd,v.tree,"AUTO-PRUNE $(h.val)", v.boxids,v.colors))
  elseif h.action == "order_prune"
    push!(vs, VizState(v.wd,v.tree,"ORDER-PRUNE $(h.val)", v.boxids,v.colors))
  elseif h.action == "orbit_prune"
    push!(vs, VizState(v.wd,v.tree,"ORBIT-PRUNE $(h.val)", v.boxids,v.colors))
  elseif h.action == "flag_skip"
    push!(vs, VizState(v.wd,v.tree,"FLAGGED FOR SKIPPING", v.boxids,v.colors))
  elseif h.action == "return"
    # for everything below h.val (a split_seq) remove coloring
    n = length(h.val)
    newcolor = Dict{Int,String}()
    for (b,clr) in v.colors
      pth=only([seq for (seq,bxid) in v.boxids if bxid==b])
      if pth == h.val
        newcolor[b]="lightblue"
      elseif common(pth, h.val) < n
        newcolor[b]=clr
      end
    end
    v = VizState(v.wd,v.tree, "Return", v.boxids, newcolor)
    push!(vs, v)
  else
    error((h.action, v))
  end
  add_history!(vs, hs[2:end])
end

function process_history(hs::Vector{History})::Vector{Pair{WiringDiagram, String}}
  @assert hs[1].action == "start"
  res, wd = Pair{WiringDiagram, String}[], WiringDiagram([],[])
  add_box!(wd, Box(view_cdict(t.coloring),[],[]))
  for h in hs[2:end]
    push!(res, add_history(wd, h))
  end
  return res
end

