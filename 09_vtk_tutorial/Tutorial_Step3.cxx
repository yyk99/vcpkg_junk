/*=========================================================================

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// First include the required header files for the VTK classes we are using.
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkConeSource.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>

int main(int, char *[]) {
    vtkNew<vtkNamedColors> colors;

    //
    // Next we create an instance of vtkConeSource and set some of its
    // properties. The instance of vtkConeSource "cone" is part of a
    // visualization pipeline (it is a source process object); it produces data
    // (output type is vtkPolyData) which other filters may process.
    //
    vtkNew<vtkConeSource> cone;
    cone->SetHeight(3.0);
    cone->SetRadius(1.0);
    cone->SetResolution(10);

    //
    // In this example we terminate the pipeline with a mapper process object.
    // (Intermediate filters such as vtkShrinkPolyData could be inserted in
    // between the source and the mapper.)  We create an instance of
    // vtkPolyDataMapper to map the polygonal data into graphics primitives. We
    // connect the output of the cone source to the input of this mapper.
    //
    vtkNew<vtkPolyDataMapper> coneMapper;
    coneMapper->SetInputConnection(cone->GetOutputPort());

    //
    // Create an actor to represent the cone. The actor orchestrates rendering
    // of the mapper's graphics primitives. An actor also refers to properties
    // via a vtkProperty instance, and includes an internal transformation
    // matrix. We set this actor's mapper to be coneMapper which we created
    // above.
    //
    vtkNew<vtkActor> coneActor;
    coneActor->SetMapper(coneMapper);
    coneActor->GetProperty()->SetColor(
        colors->GetColor3d("MistyRose").GetData());

    //
    // Create two renderers and assign actors to them. A renderer renders into
    // a viewport within the vtkRenderWindow. It is part or all of a window on
    // the screen and it is responsible for drawing the actors it has.  We also
    // set the background color here. In this example we are adding the same
    // actor to two different renderers; it is okay to add different actors to
    // different renderers as well.
    //
    vtkNew<vtkRenderer> ren1;
    ren1->AddActor(coneActor);
    ren1->SetBackground(colors->GetColor3d("RoyalBlue").GetData());

    ren1->SetViewport(0.0, 0.0, 0.5, 1.0);

    vtkNew<vtkRenderer> ren2;
    ren2->AddActor(coneActor);
    ren2->SetBackground(colors->GetColor3d("DodgerBlue").GetData());
    ren2->SetViewport(0.5, 0.0, 1.0, 1.0);

    //
    // Finally we create the render window which will show up on the screen.
    // We put our renderer into the render window using AddRenderer. We also
    // set the size to be 300 pixels by 300.
    //
    vtkNew<vtkRenderWindow> renWin;
    renWin->AddRenderer(ren1);
    renWin->AddRenderer(ren2);
    renWin->SetSize(600, 300);
    renWin->SetWindowName("Tutorial_Step3");

    //
    // Make one view 90 degrees from other.
    //
    ren1->ResetCamera();
    ren1->GetActiveCamera()->Azimuth(90);

    //
    // Now we loop over 360 degrees and render the cones each time.
    //
    for (int i = 0; i < 360; ++i) {
        // render the image
        renWin->Render();
        // rotate the active camera by one degree
        ren1->GetActiveCamera()->Azimuth(1);
        ren2->GetActiveCamera()->Azimuth(1);
    }

    return EXIT_SUCCESS;
}