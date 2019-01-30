/*
Description
    Wall function boundary condition for walls.
    Formulated from analytic solution of SA model.
    
    References:
    @verbatim
        "Modifications and Clarifications for the Implementation
         of the Spalart-Allmaras Turbulence Model"
        S.R. Allmaras,
        F.T. Johnson,
        P.R. Spalart
        Seventh International Conference on Computational Fluid Dynamics 
        (ICCFD7), 2012, pp. 1-11.
    @endverbatim
*/

#ifndef nutSAWallFunctionFvPatchScalarField_H
#define nutSAWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
           Class nutSAWallFunctionFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nutSAWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected data

        //- model constants
        scalar Bbar_;

        scalar a1_;

        scalar a2_;
        
        scalar b1_;
        
        scalar b2_;
        
        scalar c1_;
        
        scalar c2_;
        
        scalar c3_;
        
        scalar c4_;

    // Protected member functions

        //- Calculate yPLus
        virtual tmp<scalarField> calcYPlus(const scalarField& magUp) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;


public:

    //- Runtime type information
    TypeName("nutSAWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutSAWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutSAWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutSAWallFunctionFvPatchScalarField
        //  onto a new patch
        nutSAWallFunctionFvPatchScalarField
        (
            const nutSAWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutSAWallFunctionFvPatchScalarField
        (
            const nutSAWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutSAWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutSAWallFunctionFvPatchScalarField
        (
            const nutSAWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutSAWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
