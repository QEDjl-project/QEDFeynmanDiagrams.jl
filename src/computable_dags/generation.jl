function _parse_particle(name::String)
    local dir
    if startswith(name, "inc_")
        dir = Incoming()
    elseif startswith(name, "out_")
        dir = Outgoing()
    else
        throw(InvalidInputError("failed to parse particle direction from \"$name\""))
    end

    name = name[5:end]

    local species
    if startswith(name, "el")
        species = Electron()
    elseif startswith(name, "ph")
        species = Photon()
    elseif startswith(name, "po")
        species = Positron()
    else
        throw(InvalidInputError("failed to parse particle species from name \"$name\""))
    end

    name = name[4:end]

    local spin_pol
    if startswith(name, "su")
        spin_pol = SpinUp()
    elseif startswith(name, "sd")
        spin_pol = SpinDown()
    elseif startswith(name, "px")
        spin_pol = PolX()
    elseif startswith(name, "py")
        spin_pol = PolY()
    else
        throw(
            InvalidInputError(
                "failed to parse particle spin or polarization from \"$name\""
            ),
        )
    end

    name = name[4:end]

    index = parse(Int, name)
    return (dir, species, spin_pol, index)
end

function spin_or_pol(
        process::AbstractProcessDefinition,
        dir::ParticleDirection,
        species::AbstractParticleType,
        n::Int,
    )
    i = 0
    c = n
    for p in particles(process, dir)
        i += 1
        if p == species
            c -= 1
        end
        if c == 0
            break
        end
    end

    if c != 0 || n <= 0
        throw(
            InvalidInputError(
                "could not get $n-th spin/pol of $dir $species, does not exist"
            ),
        )
    end

    return spin_pols(process, dir)[i]
end

function ComputableDAGs.input_expr(
        proc::AbstractProcessDefinition, name::String, psp_symbol::Symbol
    )
    if startswith(name, "bs_")
        (dir, species, spin_pol, index) = _parse_particle(name[4:end])
        dir_str = _construction_string(dir)
        species_str = _construction_string(species)
        sp_str = _construction_string(spin_pol)

        return Meta.parse(
            "QEDFeynmanDiagrams.BaseStateInput(
                ParticleStateful($dir_str, $species_str, momentum($psp_symbol, $dir_str, $species_str, Val($index))),
                $sp_str,
            )",
        )

    elseif startswith(name, "pr_")
        index = parse(Int, name[4:end]) # get index of the virtual particle in the process

        vp = virtual_particles(proc)[index]
        return Meta.parse("QEDFeynmanDiagrams.PropagatorInput(
                              QEDFeynmanDiagrams.VirtualParticle(
                                process($psp_symbol),
                                $(_construction_string(particle_species(vp))),
                                $(vp.in_particle_contributions),
                                $(vp.out_particle_contributions)
                              ),
                              $psp_symbol
                          )")
    else
        throw(InvalidInputError("failed to parse node name \"$name\""))
    end
end

# recursion termination: base case
@inline _assemble_input_type(::Tuple{}, ::ParticleDirection) = ()

# function assembling the correct type information for the tuple of ParticleStatefuls in a phasespace point for input_type
@inline function _assemble_input_type(
        particle_types::Tuple{SPECIES_T, Vararg{AbstractParticleType}}, dir::DIR_T
    ) where {SPECIES_T <: AbstractParticleType, DIR_T <: ParticleDirection}
    return (
        AbstractParticleStateful{DIR_T, SPECIES_T},
        _assemble_input_type(particle_types[2:end], dir)...,
    )
end

function ComputableDAGs.input_type(p::AbstractProcessDefinition)
    in_t = _assemble_input_type(incoming_particles(p), Incoming())
    out_t = _assemble_input_type(outgoing_particles(p), Outgoing())
    return AbstractPhaseSpacePoint{
        typeof(p),
        <:AbstractModelDefinition,
        <:AbstractPhaseSpaceLayout,
        <:Tuple{in_t...},
        <:Tuple{out_t...},
    }
end

function Base.parse(::Type{AbstractSpinOrPolarization}, s::AbstractString)
    if s == "su"
        return SpinUp()
    end
    if s == "sd"
        return SpinDown()
    end
    if s == "px"
        return PolX()
    end
    if s == "py"
        return PolY()
    end
    throw(InvalidInputError("invalid string \"$s\" to parse to AbstractSpinOrPolarization"))
end

function _base_state_name(p::VirtualParticle)
    proc = process(p)

    dir = sum(_in_contributions(p)) == 1 ? Incoming() : Outgoing()

    # find particle in the contributions
    index = 0
    contribs = is_incoming(dir) ? _in_contributions(p) : _out_contributions(p)
    for p in particles(proc, dir)
        index += 1
        if contribs[index]
            break
        end
    end

    species = particles(proc, dir)[index]

    species_index = _species_index(proc, dir, species, index)

    spin_pol = spin_or_pol(proc, dir, species, species_index)

    return string.(
        "bs_$(_dir_str(dir))_$(_species_str(species))_",
        _spin_pol_str.(_spin_pols(spin_pol)),
        "_$(species_index)",
    )
end

# from two or three node names like "1_su-2_px"... return a single tuple of the indices and spin/pols in sorted
function _parse_node_names(name1::String, name2::String)
    split_strings_1 = split.(split(name1, "-"), "_")
    split_strings_2 = split.(split(name2, "-"), "_")

    return tuple(
        # TODO: could use merge sort since the sub lists are sorted already
        sort(
            [
                tuple.(
                    parse.(Int, getindex.(split_strings_1, 1)),
                    parse.(AbstractSpinOrPolarization, getindex.(split_strings_1, 2)),
                )...,
                tuple.(
                    parse.(Int, getindex.(split_strings_2, 1)),
                    parse.(AbstractSpinOrPolarization, getindex.(split_strings_2, 2)),
                )...,
            ]
        )...,
    )
end
function _parse_node_names(name1::String, name2::String, name3::String)
    split_strings_1 = split.(split(name1, "-"), "_")
    split_strings_2 = split.(split(name2, "-"), "_")
    split_strings_3 = split.(split(name3, "-"), "_")

    return tuple(
        # TODO: could use merge sort since the sub lists are sorted already
        sort(
            [
                tuple.(
                    parse.(Int, getindex.(split_strings_1, 1)),
                    parse.(AbstractSpinOrPolarization, getindex.(split_strings_1, 2)),
                )...,
                tuple.(
                    parse.(Int, getindex.(split_strings_2, 1)),
                    parse.(AbstractSpinOrPolarization, getindex.(split_strings_2, 2)),
                )...,
                tuple.(
                    parse.(Int, getindex.(split_strings_3, 1)),
                    parse.(AbstractSpinOrPolarization, getindex.(split_strings_3, 2)),
                )...,
            ]
        )...,
    )
end

function _make_node_name(spin_pols::Vector)
    # spin_pols is a vector of tuples Tuple{Int, AbstractSpinOrPolarization}
    node_name = ""
    first = true
    for spin_pol_tuple in spin_pols
        if !first
            node_name *= "-"
        else
            first = false
        end
        node_name *= "$(spin_pol_tuple[1])_$(_spin_pol_str(spin_pol_tuple[2]))"
    end
    return node_name
end

"""
    _is_index_valid_combination(proc::AbstractProcessDefinition, index::Tuple)

Internal function for DAG generation. Checks for a given process and a spin/pol combination whether the spin/pol combination is
part of the process, including checking for [`QEDbase.SyncedPolarization`](@extref) and [`QEDbase.SyncedSpin`](@extref).
"""
function _is_index_valid_combination(proc::AbstractProcessDefinition, index::Tuple)
    proc_spin_pols = (incoming_spin_pols(proc)..., outgoing_spin_pols(proc)...)

    # for synced spins/pols, remember the first occurrence and its definite spin/pol, then check that later ones are the same
    synced_pols = Dict{SyncedPol, AbstractDefinitePolarization}()
    synced_spins = Dict{SyncedSpin, AbstractDefiniteSpin}()

    for (i, sp) in index
        if proc_spin_pols[i] isa AllSpin || proc_spin_pols[i] isa AllPol
            # for allspin and allpol, everything is allowed
            continue
        end
        if proc_spin_pols[i] == sp
            # sp is always definite, so if they're equal the combination is always allowed
            continue
        end
        if proc_spin_pols[i] isa SyncedSpin
            if !haskey(synced_spins, proc_spin_pols[i]) # insert first occurrence
                synced_spins[proc_spin_pols[i]] = sp
                continue
            end
            if synced_spins[proc_spin_pols[i]] == sp # otherwise, check if sp is synced
                continue
            end
            # the spin is not synced
            return false
        end

        if proc_spin_pols[i] isa SyncedPol
            if !haskey(synced_pols, proc_spin_pols[i])
                synced_pols[proc_spin_pols[i]] = sp
                continue
            end
            if synced_pols[proc_spin_pols[i]] == sp
                continue
            end
            # the pol is not synced
            return false
        end

        error("encountered unknown spin or polarization type")
    end

    return true
end

"""
    graph(proc::AbstractProcessDefinition)

Generate and return a [`ComputableDAGs.DAG`](@extref), representing the computation for the squared matrix element of this scattering process, summed over spin and polarization combinations allowed by the process.
"""
function ComputableDAGs.graph(proc::PROC) where {PROC <: AbstractProcessDefinition}
    I = number_incoming_particles(proc)
    O = number_outgoing_particles(proc)
    SPECIFIC_VP = VirtualParticle{PROC, NTuple{I, Bool}, NTuple{O, Bool}}
    particles::Vector{SPECIFIC_VP} = virtual_particles(proc)                  # virtual particles that will be input to propagator tasks

    pairs = OrderedDict(particle_pairs(particles))       # pairs to generate the pair tasks
    sort!(pairs)
    triples = sort(total_particle_triples(particles))    # triples to generate the triple tasks

    g = DAG()

    # -- Base State Tasks --
    propagated_outputs = Dict{SPECIFIC_VP, Vector{Node}}()
    for dir in (Incoming(), Outgoing())
        for species in (Electron(), Positron(), Photon())
            for index in 1:number_particles(proc, dir, species)
                p = VirtualParticle(
                    proc,
                    is_outgoing(dir) ? _invert(species) : species,
                    _momentum_contribution(proc, dir, species, index)...,
                )
                for spin_pol in _spin_pols(spin_or_pol(proc, dir, species, index))
                    # gen entry nodes
                    # names are "bs_<dir>_<species>_<spin/pol>_<index>"
                    data_node_name = "bs_$(_dir_str(dir))_$(_species_str(species))_$(_spin_pol_str(spin_pol))_$(index)"

                    data_in = insert_node!(g, DataTask(0), data_node_name)

                    # generate initial base_state tasks
                    compute_base_state = insert_node!(g, ComputeTask_BaseState())

                    data_out = insert_node!(
                        g,
                        DataTask(0),
                        "$(_total_index(proc, dir, species, index))_$(_spin_pol_str(spin_pol))",
                    )

                    insert_edge!(g, data_in, compute_base_state)
                    insert_edge!(g, compute_base_state, data_out)

                    if !haskey(propagated_outputs, p)
                        propagated_outputs[p] = Vector{Node}()
                    end
                    push!(propagated_outputs[p], data_out)
                end
            end
        end
    end

    # -- Propagator Tasks --
    propagator_task_outputs = Dict()
    vp_index = 0
    for vp in virtual_particles(proc)
        vp_index += 1

        data_node_name = "pr_$vp_index"

        data_in = insert_node!(g, DataTask(0), data_node_name)
        compute_vp_propagator = insert_node!(g, ComputeTask_Propagator())
        data_out = insert_node!(g, DataTask(0))

        insert_edge!(g, data_in, compute_vp_propagator)
        insert_edge!(g, compute_vp_propagator, data_out)

        propagator_task_outputs[vp] = data_out
    end

    # -- Pair Tasks --
    for (product_particle, input_particle_vector) in pairs
        propagated_outputs[product_particle] = Vector{Node}()

        # make a dictionary of vectors to collect the outputs depending on spin/pol configs of the input particles
        N = _number_contributions(product_particle)
        pair_output_nodes_by_spin_pol = Dict{
            NTuple{N, Tuple{Int, AbstractSpinOrPolarization}}, Vector{DataTaskNode},
        }()

        for input_particles in input_particle_vector
            # input_particles is a tuple of first and second particle
            particles_data_out_nodes = (
                propagated_outputs[input_particles[1]],
                propagated_outputs[input_particles[2]],
            )

            for in_nodes in Iterators.product(particles_data_out_nodes...)
                # get the spin/pol config of the input particles from the data_out names
                index = _parse_node_names(in_nodes[1].name, in_nodes[2].name)
                # index is a tuple of tuples, containing the particle index and their definite spin/pol
                if !_is_index_valid_combination(proc, index)
                    # skip this pair creation if the spin/pol combination doesn't exist
                    continue
                end

                # make the compute pair nodes for every combination of the found input_particle_nodes to get all spin/pol combinations
                negate = relative_sign_pair(input_particles[1], input_particles[2])

                compute_pair = if negate
                    insert_node!(g, ComputeTask_PairNegated())
                else
                    insert_node!(g, ComputeTask_Pair())
                end
                pair_data_out = insert_node!(g, DataTask(0))

                insert_edge!(
                    g, in_nodes[1], compute_pair, _edge_index_from_vp(input_particles[1])
                )
                insert_edge!(
                    g, in_nodes[2], compute_pair, _edge_index_from_vp(input_particles[2])
                )
                insert_edge!(g, compute_pair, pair_data_out)

                if !haskey(pair_output_nodes_by_spin_pol, index)
                    pair_output_nodes_by_spin_pol[index] = Vector()
                end
                push!(pair_output_nodes_by_spin_pol[index], pair_data_out)
            end
        end

        propagator_node = propagator_task_outputs[product_particle]

        for (index, nodes_to_sum) in pair_output_nodes_by_spin_pol
            compute_pairs_sum = insert_node!(
                g, ComputeTask_CollectPairs()
            )

            data_pairs_sum = insert_node!(g, DataTask(0))
            compute_propagated = insert_node!(g, ComputeTask_PropagatePairs())
            # give this out node the correct name
            data_out_propagated = insert_node!(g, DataTask(0), _make_node_name([index...]))

            for node in nodes_to_sum
                insert_edge!(g, node, compute_pairs_sum, 2)
            end

            insert_edge!(g, compute_pairs_sum, data_pairs_sum)

            insert_edge!(g, propagator_node, compute_propagated, 1)
            insert_edge!(g, data_pairs_sum, compute_propagated, 2)

            insert_edge!(g, compute_propagated, data_out_propagated)

            push!(propagated_outputs[product_particle], data_out_propagated)
        end
    end

    # -- Triples --
    triples_results = Dict()
    for (ph, el, po) in triples    # for each triple (each "diagram")
        photons = propagated_outputs[ph]
        electrons = propagated_outputs[el]
        positrons = propagated_outputs[po]

        for (a, b, c) in Iterators.product(photons, electrons, positrons) # for each spin/pol config of each part
            index = _parse_node_names(a.name, b.name, c.name)
            if !_is_index_valid_combination(proc, index)
                # skip this triple creation if the spin/pol combination doesn't exist, same as for pairs
                continue
            end

            negate = relative_sign_triple(el, po, ph)

            compute_triples = if negate
                insert_node!(g, ComputeTask_TripleNegated())
            else
                insert_node!(g, ComputeTask_Triple())
            end
            data_triples = insert_node!(g, DataTask(0))

            insert_edge!(g, a, compute_triples, _edge_index_from_species(Photon())) # first argument photons
            insert_edge!(g, b, compute_triples, _edge_index_from_species(Electron())) # second argument electrons
            insert_edge!(g, c, compute_triples, _edge_index_from_species(Positron())) # third argument positrons

            insert_edge!(g, compute_triples, data_triples)

            if !haskey(triples_results, index)
                triples_results[index] = Vector{DataTaskNode}()
            end
            push!(triples_results[index], data_triples)
        end
    end

    # -- Collect Triples --
    collected_triples = Vector{DataTaskNode}()
    for (index, results) in triples_results
        compute_collect_triples = insert_node!(
            g, ComputeTask_CollectTriples()
        )
        data_collect_triples = insert_node!(g, DataTask(0))

        for triple in results
            insert_edge!(g, triple, compute_collect_triples)
        end
        insert_edge!(g, compute_collect_triples, data_collect_triples)

        push!(collected_triples, data_collect_triples)
    end

    # Finally, abs2 sum over spin/pol configurations
    compute_total_result = insert_node!(
        g, ComputeTask_SpinPolCumulation()
    )
    for finished_triple in collected_triples
        insert_edge!(g, finished_triple, compute_total_result)
    end

    final_data_out = insert_node!(g, DataTask(0))
    insert_edge!(g, compute_total_result, final_data_out)
    return g
end
