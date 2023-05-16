#/usr/bin/env bash
_uno_ampl_completions() 
{
    local cur prev opts base
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    #  The basic options to complete.
    opts="-globalization_mechanism -globalization_strategy -constraint_relaxation_strategy -subproblem -preset"

    #  Complete the arguments to some of the basic commands.
    case "${prev}" in
        -globalization_mechanism)
			local globalization_mechanisms="TR LS"
            COMPREPLY=( $(compgen -W "${globalization_mechanisms}" -- ${cur}) )
            return 0
            ;;
        -globalization_strategy)
			local globalization_strategies="l1_merit leyffer_filter_strategy waechter_filter_strategy"
            COMPREPLY=( $(compgen -W "${globalization_strategies}" -- ${cur}) )
            return 0
            ;;
        -constraint_relaxation_strategy)
			local constraint_relaxation_strategies="feasibility_restoration l1_relaxation"
            COMPREPLY=( $(compgen -W "${constraint_relaxation_strategies}" -- ${cur}) )
            return 0
            ;;
		-subproblem)
			local subproblems="QP LP primal_dual_interior_point"
            COMPREPLY=( $(compgen -W "${subproblems}" -- ${cur}) )
            return 0
            ;;
        -preset)
			local presets="filtersqp ipopt byrd"
            COMPREPLY=( $(compgen -W "${presets}" -- ${cur}) )
            return 0
            ;;
        *)
        ;;
    esac

	if [[ ${cur} == -* ]] ; then
		COMPREPLY=($(compgen -W "${opts}" -- ${cur}))
		return 0
	else
		_filedir
    fi
}
complete -F _uno_ampl_completions uno_ampl
