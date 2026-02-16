const PROJECT_ROOT = dirname(dirname(@__FILE__))
cd(PROJECT_ROOT)

using Revise
using EnzymeRates, Test, Random, LinearAlgebra
using ModelContextProtocol

const _EVAL_FILE = joinpath(PROJECT_ROOT, ".mcp_eval.jl")

function _trim_stacktrace(bt_str::String)
    lines = split(bt_str, '\n')
    # Keep frames until we hit ModelContextProtocol internals
    cutoff = findfirst(l -> contains(l, "ModelContextProtocol"), lines)
    if cutoff !== nothing && cutoff > 1
        return join(lines[1:cutoff-1], '\n')
    end
    return bt_str
end

function _capture_eval(code::String)
    out_pipe = Pipe()
    err_pipe = Pipe()
    result = nothing
    err = nothing
    redirect_stdio(; stdout=out_pipe, stderr=err_pipe) do
        try
            # Write to a temp file in PROJECT_ROOT so that include()
            # inside the code resolves relative paths correctly.
            write(_EVAL_FILE, code)
            result = cd(PROJECT_ROOT) do
                Base.include(Main, _EVAL_FILE)
            end
        catch e
            err = (e, catch_backtrace())
        finally
            rm(_EVAL_FILE; force=true)
        end
    end
    close(out_pipe.in)
    close(err_pipe.in)
    out_str = read(out_pipe, String)
    err_str = read(err_pipe, String)
    if err !== nothing
        buf = IOBuffer()
        print(buf, "ERROR: ")
        showerror(buf, err[1], err[2])
        return CallToolResult(
            content=[Dict("type" => "text", "text" => _trim_stacktrace(String(take!(buf))))],
            is_error=true,
        )
    end
    parts = String[]
    isempty(out_str) || push!(parts, out_str)
    isempty(err_str) || push!(parts, err_str)
    result_repr = repr(result)
    if result_repr != "nothing"
        push!(parts, result_repr)
    end
    text = isempty(parts) ? "nothing" : join(parts, "\n")
    return TextContent(text=text)
end

server = mcp_server(
    name="julia-repl",
    version="0.1.0",
    description="Julia REPL for EnzymeRates development",
    tools=[
        MCPTool(
            name="exec_julia",
            description="Execute Julia code in a persistent session with EnzymeRates loaded. Revise.jl is active so source edits are picked up automatically.",
            parameters=[
                ToolParameter(
                    name="code",
                    description="Julia code to evaluate",
                    type="string",
                    required=true,
                ),
            ],
            handler=function (params)
                code = get(params, "code", "")
                return _capture_eval(code)
            end,
        ),
    ],
)

start!(server)
