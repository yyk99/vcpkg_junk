#include <vtkActor.h>
#include <vtkAxesActor.h>
#include <vtkCamera.h>
#include <vtkCaptionActor2D.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkTextProperty.h>
#include <vtkTransform.h>

int main(int, char *[]) {
    vtkNew<vtkNamedColors> colors;

    vtkNew<vtkSphereSource> sphereSource;
    sphereSource->SetCenter(0.0, 0.0, 0.0);
    sphereSource->SetRadius(0.5);

    // create a mapper
    vtkNew<vtkPolyDataMapper> sphereMapper;
    sphereMapper->SetInputConnection(sphereSource->GetOutputPort());

    // create an actor
    vtkNew<vtkActor> sphereActor;
    sphereActor->SetMapper(sphereMapper);

    // a renderer and render window
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->SetWindowName("Axes");
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(300, 300);

    // an interactor
    vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // add the actors to the scene
    renderer->AddActor(sphereActor);
    renderer->SetBackground(colors->GetColor3d("SlateGray").GetData());

    vtkNew<vtkTransform> transform;
    transform->Translate(1.0, 0.0, 0.0);

    vtkNew<vtkAxesActor> axes;

    // The axes are positioned with a user transform
    axes->SetUserTransform(transform);

    // properties of the axes labels can be set as follows
    // this sets the x axis label to red
    // axes->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetColor(
    //  colors->GetColor3d("Red").GetData());

    // the actual text of the axis label can be changed:
    // axes->SetXAxisLabelText("test");

    renderer->AddActor(axes);

    renderer->GetActiveCamera()->Azimuth(50);
    renderer->GetActiveCamera()->Elevation(-30);

    renderer->ResetCamera();
    renderWindow->SetWindowName("Axes");
    renderWindow->Render();

    // begin mouse interaction
    renderWindowInteractor->Start();

    return EXIT_SUCCESS;
}