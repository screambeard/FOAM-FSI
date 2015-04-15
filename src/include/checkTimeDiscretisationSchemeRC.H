const IOdictionary& fvSchemes = mesh.lookupObject<IOdictionary>("fvSchemes");
const dictionary& ddtSchemes = fvSchemes.subDict("ddtSchemes");

word ddtScheme;
if(ddtSchemes.found("ddt(U)"))
{
	ddtScheme = word(ddtSchemes.lookup("ddt(U)"));
}
else
{
	ddtScheme = word(ddtSchemes.lookup("default"));	
}

if(ddtScheme!=word("EulerRC") && ddtScheme!=word("backwardRC") && ddtScheme!=word("myBackward") && ddtScheme!=word("myEuler"))
{
	FatalErrorIn("icoDySymFoamP checkTimeDiscretisationScheme: ")
        << " ddt(U) scheme is currently set to " << ddtScheme
        << ". This solver only works with ddt(U) scheme = EulerRC or backwardRC."
        << abort(FatalError);
}
