
/* *********************************************************   
 *          
 *   Last modified by:    $Author: rrossi $
 *   Date:                $Date: 2003/10/06 12:53:21 $
 *   Revision:            $Revision: 1.15 $
 *
 * ***********************************************************/


/* System includes */
#include <set>


/* External includes */


/* Project includes */
#include "../includes/define.h"
#include "../includes/dof.h"
#include "../includes/model.h"
#include "../includes/scheme.h"
#include "../includes/update.h"

namespace Kratos
{
  
Update::Update(Model::Pointer pNewModel): mpModel(pNewModel)
{
}

Update::~Update()
{
}

//*******************************************************************************
//* - Update of all the relevant variables and derivatives, done using directly the 
//*   solution vector. Elements are NOT used in this operation.
//********************************************************************************
void Update::UpdateVariables(Vector<double>& rX,
						     std::set<Dof::Pointer,ComparePDof>& rDofSet)
{
}






} // Namespace Kratos
