#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/polygon/voronoi.hpp>

namespace bg = boost::geometry;
typedef bg::model::d2::point_xy<double> point2d_t;
typedef bg::model::polygon<point2d_t> polygon2d_t;

int main(int, char *[]) {

    vtkNew<vtkNamedColors> colors;

    // Set the background color.
    std::array<unsigned char, 4> bkg{{26, 51, 102, 255}};
    colors->SetColor("BkgColor", bkg.data());

    polygon2d_t wkt_polygon;
    bg::read_wkt("POLYGON((0 0, 0 2, 2 2, 2 0, 0 0))", wkt_polygon);

    // vtkPoints represents 3D points. The data model for vtkPoints is an array
    // of vx-vy-vz triplets accessible by (point or cell) id.
    vtkNew<vtkPoints> points;
    points->SetNumberOfPoints(wkt_polygon.outer().size());

    // vtkCellArray is a supporting object that explicitly represents cell
    // connectivity.
    // The cell array structure is a raw integer list of the form:
    // (n,id1,id2,...,idn, n,id1,id2,...,idn, ...) where n is the number of
    // points in the cell, and id is a zero-offset index into an associated
    // point list.
    vtkNew<vtkCellArray> lines;
    lines->InsertNextCell(wkt_polygon.outer().size());
    for (int i = 0; i != wkt_polygon.outer().size() - 1; ++i) {
        auto &v = wkt_polygon.outer()[i];
        points->SetPoint(i, v.get<0>(), v.get<1>(), 0.0);
        lines->InsertCellPoint(i);
    }
    lines->InsertCellPoint(0);

    // vtkPolyData is a data object that is a concrete implementation of
    // vtkDataSet.
    // vtkPolyData represents a geometric structure consisting of vertices,
    // lines, polygons, and/or triangle strips
    vtkNew<vtkPolyData> polygon;
    polygon->SetPoints(points);
    polygon->SetLines(lines);

    // vtkPolyDataMapper is a class that maps polygonal data (i.e., vtkPolyData)
    // to graphics primitives
    vtkNew<vtkPolyDataMapper> polygonMapper;
    polygonMapper->SetInputData(polygon);
    polygonMapper->Update();

    // Create an actor to represent the polygon. The actor orchestrates
    // rendering of the mapper's graphics primitives. An actor also refers to
    // properties via a vtkProperty instance, and includes an internal
    // transformation matrix. We set this actor's mapper to be polygonMapper
    // which we created above.
    vtkNew<vtkActor> polygonActor;
    polygonActor->SetMapper(polygonMapper);
    polygonActor->GetProperty()->SetColor(
        colors->GetColor3d("AliceBlue").GetData());

    // Create the Renderer and assign actors to it. A renderer is like a
    // viewport. It is part or all of a window on the screen and it is
    // responsible for drawing the actors it has.  We also set the
    // background color here.
    vtkNew<vtkRenderer> ren;
    ren->AddActor(polygonActor);
    ren->SetBackground(colors->GetColor3d("BkgColor").GetData());

    // Automatically set up the camera based on the visible actors.
    // The camera will reposition itself to view the center point of the actors,
    // and move along its initial view plane normal
    // (i.e., vector defined from camera position to focal point) so that all of
    // the
    // actors can be seen.
    ren->ResetCamera();

    // Finally we create the render window which will show up on the screen
    // We put our renderer into the render window using AddRenderer. We
    // also set the size to be 300 pixels by 300.
    vtkNew<vtkRenderWindow> renWin;
    renWin->SetWindowName("OrderedPolyLine");
    renWin->AddRenderer(ren);
    renWin->SetSize(300, 300);

    // The vtkRenderWindowInteractor class watches for events (e.g., keypress,
    // mouse) in the vtkRenderWindow. These events are translated into
    // event invocations that VTK understands (see VTK/Common/vtkCommand.h
    // for all events that VTK processes). Then observers of these VTK
    // events can process them as appropriate.
    vtkNew<vtkRenderWindowInteractor> iren;
    iren->SetRenderWindow(renWin);
    renWin->Render();
    iren->Initialize();
    iren->Start();

    return EXIT_SUCCESS;
}