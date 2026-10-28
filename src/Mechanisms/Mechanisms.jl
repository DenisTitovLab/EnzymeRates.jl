# ABOUTME: The Mechanisms module: reactions and mechanisms as values (raw input, canonical
# ABOUTME: graph, Sig), their naming rules, and the DSL's unchecked reader and printer.
module Mechanisms

import Unicode
using SHA: sha256

const API = (:AbstractEnzymeMechanism, :AllostericEnzymeMechanism, :AllostericMechanism,
    :AllostericRegulator, :CompetitiveInhibitor, :Counts, :Edge, :EnzymeMechanism,
    :EnzymeReaction, :Form, :Mechanism, :Met, :Metabolites, :RawAllostericMechanism,
    :RawForm, :RawMechanism, :RawStep, :RegulatorySite, :_AnyMech, :_Column, :_bound_entry,
    :_build, :_canonicalize, :_check_aligned, :_check_uniquely_decodable,
    :_chemistry_groups, :_concrete, :_edge_constants, :_entry, :_entry_names, :_excluded_by,
    :_excludes, :_exclusion_components, :_flags, :_flat_edges, :_form_name, :_form_text,
    :_group_members, :_groups_with, :_has_inactive_columns, :_has_unbalanced_cycle,
    :_in_active_law, :_is_binding, :_is_chemistry, :_is_iso, :_kept_members, :_kreg_params,
    :_metabolite_names, :_modifier_mask, :_name_index, :_names_mask, :_number_groups,
    :_ordered_columns, :_parse_reaction_block, :_raw, :_raw_form, :_raw_mechanism,
    :_raw_sites, :_reaction_expr, :_reg_K, :_runs_toward_products, :_segment_walk, :_sides,
    :_solve_tags, :_step_text, :_with, :allo_states, :allowed_catalytic_multiplicities,
    :atoms, :bound, :cat_allo_states, :catalytic, :catalytic_multiplicity,
    :compile_mechanism, :conformation, :conformation_steps, :copies, :ends, :entries,
    :given_off, :group, :has_turnover_cycle, :is_rapid_equilibrium, :ligands, :mechanism_id,
    :members, :metabolite, :metabolite_names, :metabolites, :multiplicity, :name, :ngroups,
    :present_species, :products, :re_segments, :reactants, :reg_type, :regulator,
    :regulators, :regulatory_sites, :residual, :segment_graph, :shared_catalytic_site,
    :species, :step_kind, :steps, :substrates, :tag, :taken_up)

include("types.jl")
include("accessors.jl")
include("names.jl")
include("dsl.jl")

end
