
#include "nutSAWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutSAWallFunctionFvPatchScalarField::calcYPlus
(
    const scalarField& magUp
) const
{
    const label patchi = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    const scalarField& y = rasModel.y()[patchi];
    const fvPatchScalarField& nuw = rasModel.nu().boundaryField()[patchi];

    tmp<scalarField> tyPlus(new scalarField(patch().size(), 0.0));
    scalarField& yPlus = tyPlus();

    forAll(yPlus, faceI)
    {
        scalar Re = magUp[faceI]*y[faceI]/nuw[faceI];

        scalar yp = yPlusLam_;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yp;
            scalar f = Bbar_
                     + c1_*log(sqr(yp + a1_) + sqr(b1_))
                     - c2_*log(sqr(yp + a2_) + sqr(b2_))
                     - c3_*atan2(b1_, yp+a1_)
                     - c4_*atan2(b2_, yp+a2_);
            
            yp = (Re + kappa_*yp)/(kappa_ + f);

        } while (mag(yp - yPlusLast) > 1e-8 && ++iter < 1000);

        yPlus[faceI] = max(0.0, yp);
    }

    return tyPlus;
}


tmp<scalarField>
nutSAWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];
    const scalarField magUp = mag(Uw.patchInternalField() - Uw);
    const scalarField& nuw = rasModel.nu().boundaryField()[patchI];

    tmp<scalarField> tyPlus = calcYPlus(magUp);
    scalarField& yPlus = tyPlus();

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw();

    forAll(yPlus, faceI)
    {
        scalar f = Bbar_
                 + c1_*log(sqr(yPlus[faceI] + a1_) + sqr(b1_))
                 - c2_*log(sqr(yPlus[faceI] + a2_) + sqr(b2_))
                 - c3_*atan2(b1_, yPlus[faceI]+a1_)
                 - c4_*atan2(b2_, yPlus[faceI]+a2_);

        nutw[faceI] = max(nuw[faceI]*(yPlus[faceI]/f - 1.0), scalar(0.0));
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutSAWallFunctionFvPatchScalarField::nutSAWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    Bbar_(5.0333908790505579),
    a1_(8.148221580024245),
    a2_(-6.9287093849022945),
    b1_(7.4600876082527945),
    b2_(7.468145790401841),
    c1_(2.5496773539754747),
    c2_(1.3301651588535228),
    c3_(3.599459109332379),
    c4_(3.6397531868684494)
{}


nutSAWallFunctionFvPatchScalarField::nutSAWallFunctionFvPatchScalarField
(
    const nutSAWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Bbar_(ptf.Bbar_),
    a1_(ptf.a1_),
    a2_(ptf.a2_),
    b1_(ptf.b1_),
    b2_(ptf.b2_),
    c1_(ptf.c1_),
    c2_(ptf.c2_),
    c3_(ptf.c3_),
    c4_(ptf.c4_)
{}


nutSAWallFunctionFvPatchScalarField::nutSAWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    Bbar_(dict.lookupOrDefault<scalar>("Bbar", 5.0333908790505579)),
    a1_(dict.lookupOrDefault<scalar>("a1", 8.148221580024245)),
    a2_(dict.lookupOrDefault<scalar>("a2", -6.9287093849022945)),
    b1_(dict.lookupOrDefault<scalar>("b1", 7.4600876082527945)),
    b2_(dict.lookupOrDefault<scalar>("b2", 7.468145790401841)),
    c1_(dict.lookupOrDefault<scalar>("c1", 2.5496773539754747)),
    c2_(dict.lookupOrDefault<scalar>("c2", 1.3301651588535228)),
    c3_(dict.lookupOrDefault<scalar>("c3", 3.599459109332379)),
    c4_(dict.lookupOrDefault<scalar>("c4", 3.6397531868684494))
{}


nutSAWallFunctionFvPatchScalarField::nutSAWallFunctionFvPatchScalarField
(
    const nutSAWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf),
    Bbar_(wfpsf.Bbar_),
    a1_(wfpsf.a1_),
    a2_(wfpsf.a2_),
    b1_(wfpsf.b1_),
    b2_(wfpsf.b2_),
    c1_(wfpsf.c1_),
    c2_(wfpsf.c2_),
    c3_(wfpsf.c3_),
    c4_(wfpsf.c4_)
{}


nutSAWallFunctionFvPatchScalarField::nutSAWallFunctionFvPatchScalarField
(
    const nutSAWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    Bbar_(wfpsf.Bbar_),
    a1_(wfpsf.a1_),
    a2_(wfpsf.a2_),
    b1_(wfpsf.b1_),
    b2_(wfpsf.b2_),
    c1_(wfpsf.c1_),
    c2_(wfpsf.c2_),
    c3_(wfpsf.c3_),
    c4_(wfpsf.c4_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
nutSAWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];
    const scalarField magUp = mag(Uw.patchInternalField() - Uw);

    return calcYPlus(magUp);
}


void nutSAWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    os.writeKeyword("Bbar") << Bbar_ << token::END_STATEMENT << nl;
    os.writeKeyword("a1") << a1_ << token::END_STATEMENT << nl;
    os.writeKeyword("a2") << a2_ << token::END_STATEMENT << nl;
    os.writeKeyword("b1") << b1_ << token::END_STATEMENT << nl;
    os.writeKeyword("b2") << b2_ << token::END_STATEMENT << nl;
    os.writeKeyword("c1") << c1_ << token::END_STATEMENT << nl;
    os.writeKeyword("c2") << c2_ << token::END_STATEMENT << nl;
    os.writeKeyword("c3") << c3_ << token::END_STATEMENT << nl;
    os.writeKeyword("c4") << c4_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nutSAWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
