mutable struct CustomModel <: Model
  params
  pdflist::Array{NLogPDF}
  nll::NLogLikelihood
  function CustomModel()
    new( [], 
        Array{NLogPDF}(undef, 0),
        NLogLikelihood() )
  end
end

function add_parameters!( m::CustomModel, params )
  for p in params
    push!(m.params, p)
  end
end

function add_nlogpdfs!( m::CustomModel, pdflist::Array{NLogPDF} )
  for pdf in pdflist
    push!(m.pdflist, pdf)
  end
end

function minimize!( m::CustomModel )
  likelihood = NLogLikelihood( m.pdflist )
  add_likelihood!(m, likelihood )
  optimize_model!( m, likelihood )
end
