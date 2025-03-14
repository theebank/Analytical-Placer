/*
 * Copyright 2019-2023 University of Toronto
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Authors: Mario Badr, Sameh Attia, Tanner Young-Schultz and Vaughn Betz
 */

/**
 * @file
 *
 * This example shows you how to create an application using the EZGL library.
 */

#include <chrono>
#include <thread>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <utility>
#include <map>
#include <limits>
#include <algorithm>
#include "ezgl_code_docs_example/ezgl_example/ezgl/application.hpp"
#include "ezgl_code_docs_example/ezgl_example/ezgl/graphics.hpp"
#include "analyticalPlacer.hpp"

using namespace std;
// FUNCTION DECLARATIONS

/**
 * Draw to the main canvas using the provided graphics object. Runs every time graphics are refreshed/image zooms or pans
 *
 * The graphics object expects that x and y values will be in the main canvas' world coordinate system.
 */
void draw_main_canvas(ezgl::renderer *g);

/**
 * Initial Setup is a mandatory function for any EZGL application, and is run whenever a window is opened.
 */
void initial_setup(ezgl::application *application, bool new_window);

/**
 * DRAWING HELPER FUNCTIONS
 *
 * draw_main_canvas helper functions. Example functions that draw different things.
 */
void draw_rectangle_example(ezgl::renderer *g);
void draw_arc_example(ezgl::renderer *g);
void rotated_text_example(ezgl::renderer *g);
void draw_poly_example(ezgl::renderer *g);
void draw_text_example(ezgl::renderer *g);
void draw_line_example(ezgl::renderer *g);
void screen_coordinates_example(ezgl::renderer *g);
void draw_png_example(ezgl::renderer *g);

void draw_grid(ezgl::renderer *g);
void draw_blocks(ezgl::renderer *g);
void draw_connections(ezgl::renderer *g);
void draw_HPWL(ezgl::renderer *g);
void draw_anchors(ezgl::renderer *g);

/**
 * UI CALLBACK FUNCTIONS
 *
 * These are example callback functions for the UI elements
 */
void animate_button_cbk(GtkWidget *widget, ezgl::application *application);
void test_button_cbk(GtkWidget *widget, ezgl::application *application);
void combo_box_cbk(GtkComboBoxText *self, ezgl::application *app);
void delete_combo_box_cbk(GtkWidget *widget, ezgl::application *application);
void create_dialog_button_cbk(GtkWidget *widget, ezgl::application *application);
void create_mssg_button_cbk(GtkWidget *widget, ezgl::application *application);
void dialog_cbk(GtkDialog *self, gint response_id, ezgl::application *app);

/**
 * EVENT CALLBACK FUNCTIONS
 *
 * These functions run whenever their corresponding event (key press, mouse move, or mouse click) occurs.
 */
void act_on_mouse_press(ezgl::application *application, GdkEventButton *event, double x, double y);
void act_on_mouse_move(ezgl::application *application, GdkEventButton *event, double x, double y);
void act_on_key_press(ezgl::application *application, GdkEventKey *event, char *key_name);

static ezgl::rectangle initial_world{{0, 0}, 1100, 1150};

analyticalPlacer ap;
int CellDisplacement = 0;
bool anchors;

/**
 * The start point of the program.
 *
 * This function initializes an ezgl application and runs it.
 *
 * @param argc The number of arguments provided.
 * @param argv The arguments as an array of c-strings.
 *
 * @return the exit status of the application run.
 */

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        std::cerr << "usage: " << argv[0] << " <test circuit Num> <weight factor between fixed node and movables> <0 for Parts 1/2, 1 for Part 3.1, 2 for Part 3.2>\n";
        return 1;
    }

    double weightfactor = std::atof(argv[2]);
    ezgl::application::settings settings;

    // Path to the "main.ui" file that contains an XML description of the UI.
    // Edit this file with glade if you want to change the UI layout
    settings.main_ui_resource = "main.ui";

    // Note: the "main.ui" file has a GtkWindow called "MainWindow".
    settings.window_identifier = "MainWindow";

    // Note: the "main.ui" file has a GtkDrawingArea called "MainCanvas".
    settings.canvas_identifier = "MainCanvas";

    // Create our EZGL application.
    ezgl::application application(settings);

    // Set some parameters for the main sub-window (MainCanvas), where
    // visualization graphics are draw. Set the callback function that will be
    // called when the main window needs redrawing, and define the (world)
    // coordinate system we want to draw in.
    application.add_canvas("MainCanvas", draw_main_canvas, initial_world);

    string filenum = argv[1];
    ap.analyticalPlacementv2("./tests/cct" + filenum + ".txt", weightfactor);
    int number;
    try
    {
        number = stoi(argv[3]);
    }
    catch (const invalid_argument &e)
    {
        cerr << "Cannot convert to int" << endl;
    }
    if (number == 1)
    {
        CellDisplacement = ap.daravspreading();
        ap.HPWL = ap.calculateHPWL();
    }
    else if (number == 2)
    {
        CellDisplacement = ap.daravspreading();
        ap.analyticalPlacementAnchored("./tests/cct" + filenum + ".txt", weightfactor);
        anchors = true;
    }

    // Run the application until the user quits.
    // This hands over all control to the GTK runtime---after this point
    // you will only regain control based on callbacks you have setup.
    // Three callbacks can be provided to handle mouse button presses,
    // mouse movement and keyboard button presses in the graphics area,
    // respectively. Also, an initial_setup function can be passed that will
    // be called before the activation of the application and can be used
    // to create additional buttons, initialize the status message, or
    // connect added widgets to their callback functions.
    // Those callbacks are optional, so we can pass nullptr if
    // we don't need to take any action on those events
    return application.run(initial_setup, act_on_mouse_press, act_on_mouse_move, act_on_key_press);
}

/**
 * The redrawing function for still pictures
 */
void draw_main_canvas(ezgl::renderer *g)
{
    draw_grid(g);
    draw_connections(g);
    draw_blocks(g);
    draw_HPWL(g);
    if (anchors)
    {
        draw_anchors(g);
    }
}

/**
 * Function called before the activation of the application
 * Can be used to create additional buttons, initialize the status message,
 * or connect added widgets to their callback functions
 */
void initial_setup(ezgl::application *application, bool /*new_window*/)
{
    // Update the status bar message
    application->update_message("EZGL Application");

    // Setting our starting row for insertion at 6 (Default zoom/pan buttons created by EZGL take up first five rows);
    // We will increment row each time we insert a new element.
    int row = 6;

    application->create_label(row++, "EXAMPLE UI ELEMENTS: ");

    // Create the Animate button and link it with animate_button callback fn.
    application->create_button("Animate", row++, animate_button_cbk);

    // Create a Test button and link it with test_button callback fn.
    application->create_button("Test", row++, test_button_cbk);

    application->create_label(row++, "Test Combo Box:");

    // Creating example combo box with options Yes, No, Maybe, connected to combo_box_cbk
    application->create_combo_box_text(
        "TestComboBox",
        row++,
        combo_box_cbk,
        {"YES", "NO", "MAYBE"});

    application->create_button("Delete Combo Box", row++, delete_combo_box_cbk);

    application->create_button("Create Dialog", row++, create_dialog_button_cbk);

    application->create_button("Create Popup Mssg", row++, create_mssg_button_cbk);
}

/**
 * Draw some rectangles with different colors
 */
void draw_rectangle_example(ezgl::renderer *g)
{
    const float rectangle_width = 50;
    const float rectangle_height = rectangle_width;
    const ezgl::point2d start_point(150, 30);
    ezgl::rectangle color_rectangle = {start_point, rectangle_width, rectangle_height};

    // Some of the available colors, a complete list is in ezgl/include/color.hpp
    ezgl::color color_indicies[] = {
        ezgl::GREY_55,
        ezgl::GREY_75,
        ezgl::WHITE,
        ezgl::BLACK,
        ezgl::BLUE,
        ezgl::GREEN,
        ezgl::YELLOW,
        ezgl::CYAN,
        ezgl::RED,
        ezgl::DARK_GREEN,
        ezgl::MAGENTA};

    // format text font and color
    g->set_color(ezgl::BLACK);
    g->format_font("monospace", ezgl::font_slant::normal, ezgl::font_weight::normal, 10);

    // draw text
    g->draw_text({110, color_rectangle.center_y()}, "colors", 2 * (start_point.x - 110), rectangle_height);

    for (size_t i = 0; i < sizeof(color_indicies) / sizeof(color_indicies[0]); ++i)
    {
        // Change the color of next draw calls
        g->set_color(color_indicies[i]);

        // Draw filled in rectangles
        g->fill_rectangle(color_rectangle);

        // Increment the start point
        color_rectangle += {rectangle_width, 0};
    }

    // Draw text
    g->draw_text({400, color_rectangle.center_y()}, "fill_rectangle", DBL_MAX, rectangle_height);

    /* Draw some rectangles with RGB triplet colors and alpha (transparency) */

    // Hack to make the colors change once per second
    std::srand(time(0));

    for (size_t i = 0; i < 3; ++i)
    {
        // Increment the start point
        color_rectangle += {rectangle_width, 0};

        // Change the next draw calls color. rgb and alpha values range from 0 to 255
        g->set_color(std::rand() % 256, std::rand() % 256, std::rand() % 256, 255);

        // Draw filled in rectangles
        g->fill_rectangle(color_rectangle);
    }

    /* Draw a black border rectangle */

    // Change the next draw calls color to black
    g->set_color(ezgl::BLACK);

    // Change the next draw calls line width
    g->set_line_width(1);

    // Draw a rectangle bordering all the drawn rectangles
    g->draw_rectangle(start_point, color_rectangle.top_right());
}

/**
 * Draw some example lines, shapes, and arcs
 */
void draw_arc_example(ezgl::renderer *g)
{
    float radius = 50;

    // Draw solid line
    g->set_color(ezgl::BLACK);
    g->draw_text({250, 150}, "draw_line", 150.0, DBL_MAX);
    g->set_line_dash(ezgl::line_dash::none);
    g->draw_line({200, 120}, {200, 200});

    // Draw dashed line
    g->set_line_dash(ezgl::line_dash::asymmetric_5_3);
    g->draw_line({300, 120}, {300, 200});

    // Draw elliptic arc
    g->set_color(ezgl::MAGENTA);
    g->draw_text({450, 160}, "draw_elliptic_arc", 150.0, DBL_MAX);
    g->draw_elliptic_arc({550, 160}, 30, 60, 90, 270);

    // Draw filled in elliptic arc
    g->draw_text({700, 160}, "fill_elliptic_arc", 150.0, DBL_MAX);
    g->fill_elliptic_arc({800, 160}, 30, 60, 90, 270);

    // Draw arcs
    g->set_color(ezgl::BLUE);
    g->draw_text({190, 300}, "draw_arc", radius * 2, 150);
    g->draw_arc({190, 300}, radius, 0, 270);
    g->draw_arc({300, 300}, radius, 0, -180);

    // Draw filled in arcs
    g->fill_arc({410, 300}, radius, 90, -90);
    g->fill_arc({520, 300}, radius, 0, 360);
    g->set_color(ezgl::BLACK);
    g->draw_text({520, 300}, "fill_arc", radius * 2, 150);
    g->set_color(ezgl::BLUE);
    g->fill_arc({630, 300}, radius, 90, 180);
    g->fill_arc({740, 300}, radius, 90, 270);
    g->fill_arc({850, 300}, radius, 90, 30);
}

/**
 * Draw some rotated text
 */
void rotated_text_example(ezgl::renderer *g)
{
    const float textsquare_width = 200;

    ezgl::rectangle textsquare = {{100, 400}, textsquare_width, textsquare_width};

    g->set_color(ezgl::BLUE);
    g->draw_rectangle(textsquare);

    g->set_color(ezgl::GREEN);
    g->draw_rectangle(textsquare.center(), {textsquare.right(), textsquare.top()});
    g->draw_rectangle({textsquare.left(), textsquare.bottom()}, textsquare.center());

    g->set_color(ezgl::RED);
    g->draw_line({textsquare.left(), textsquare.bottom()}, {textsquare.right(), textsquare.top()});
    g->draw_line({textsquare.left(), textsquare.top()}, {textsquare.right(), textsquare.bottom()});

    g->set_color(0, 0, 0, 100);
    g->set_font_size(14);
    g->draw_text({textsquare.center_x(), textsquare.bottom()}, "0 degrees", textsquare.width(), textsquare.height());

    g->set_text_rotation(90);
    g->draw_text({textsquare.right(), textsquare.center_y()}, "90 degrees", textsquare.width(), textsquare.height());

    g->set_text_rotation(180);
    g->draw_text({textsquare.center_x(), textsquare.top()}, "180 degrees", textsquare.width(), textsquare.height());

    g->set_text_rotation(270);
    g->draw_text({textsquare.left(), textsquare.center_y()}, "270 degrees", textsquare.width(), textsquare.height());

    g->set_text_rotation(45);
    g->draw_text(textsquare.center(), "45 degrees", textsquare.width(), textsquare.height());

    g->set_text_rotation(135);
    g->draw_text(textsquare.center(), "135 degrees", textsquare.width(), textsquare.height());

    // It is probably a good idea to set text rotation back to zero,
    g->set_text_rotation(0);
}

/**
 * Draw some Polygons
 */
void draw_poly_example(ezgl::renderer *g)
{
    g->set_font_size(10);
    g->set_color(ezgl::RED);

    // Draw a triangle
    g->fill_poly({{500, 400}, {440, 480}, {560, 480}});

    // Draw a 4-point polygon
    g->fill_poly({{700, 400}, {650, 480}, {750, 480}, {800, 400}});

    g->set_color(ezgl::BLACK);
    g->draw_text({500, 450}, "fill_poly", 80.0, DBL_MAX);
    g->draw_text({725, 440}, "fill_poly", 100.0, DBL_MAX);

    g->set_color(ezgl::DARK_GREEN);
    g->set_line_dash(ezgl::line_dash::none);
    ezgl::rectangle rect = {{350, 550}, {650, 670}};
    g->draw_text(rect.center(), "draw_rectangle", rect.width(), rect.height());
    g->draw_rectangle(rect);

    /* Draw some semi-transparent primitives */
    g->set_font_size(10);

    g->set_color(255, 0, 0, 255);
    g->fill_rectangle({1000, 400}, {1050, 800});

    g->set_color(0, 0, 255, 255);
    g->fill_rectangle({1000 + 50, 400}, {1050 + 50, 800});

    g->set_color(0, 255, 0, 255 / 2); // 50% transparent
    g->fill_rectangle({1000 + 25, 400 - 100}, {1050 + 25, 800 - 200});

    g->set_color(255, 100, 255, 255 / 2);
    g->fill_poly({{465, 380}, {400, 450}, {765, 450}, {850, 380}});

    g->set_color(100, 100, 255, 255 / 3);
    g->fill_poly({{550, 420}, {475, 500}, {875, 500}});

    g->set_color(ezgl::BLACK);
    g->set_text_rotation(90);
    g->draw_text({1000 - 50, 500}, "Partially transparent polys", 500, DBL_MAX);
    g->set_text_rotation(0);
}

/**
 * Draw some example text, with the bounding box functions
 */
void draw_text_example(ezgl::renderer *g)
{

    const float text_example_width = 800;
    const int num_lines = 2;
    const int max_strings_per_line = 3;
    const int num_strings_per_line[num_lines] = {3, 2};

    const char *const line_text[num_lines][max_strings_per_line] = {
        {"8 Point Text",
         "12 Point Text",
         "18 Point Text"},
        {"24 Point Text",
         "32 Point Text"}};

    const int text_sizes[num_lines][max_strings_per_line] = {
        {8, 12, 15},
        {24, 32}};

    g->set_color(ezgl::BLACK);
    g->set_line_dash(ezgl::line_dash::asymmetric_5_3);

    for (int i = 0; i < num_lines; ++i)
    {
        ezgl::rectangle text_bbox = {{100., 710. + i * 60.}, text_example_width / num_strings_per_line[i], 60.};

        for (int j = 0; j < num_strings_per_line[i]; ++j)
        {
            g->set_font_size(text_sizes[i][j]);
            g->draw_text(text_bbox.center(), line_text[i][j], text_bbox.width(), text_bbox.height());
            g->draw_rectangle(text_bbox);
            text_bbox = {{text_bbox.left() + text_example_width / num_strings_per_line[i], text_bbox.bottom()}, text_bbox.width(), text_bbox.height()};
        }
    }
}

/**
 * Draw wide lines with different end shapes
 */
void draw_line_example(ezgl::renderer *g)
{
    g->set_font_size(10);

    for (int i = 0; i <= 2; ++i)
    {
        double offsetY = 50 * i;

        g->set_horiz_justification(ezgl::justification::left);

        if (i == 0)
        {
            g->set_color(ezgl::BLACK);
            g->set_line_cap(ezgl::line_cap::butt);   // Butt ends
            g->set_line_dash(ezgl::line_dash::none); // Solid line
            g->draw_text({950, 920 + offsetY}, "Butt ends, opaque", 400, DBL_MAX);
        }

        else if (i == 1)
        {
            g->set_color(ezgl::GREEN, 255 * 2 / 3);  // Green line that is 33% transparent)
            g->set_line_cap(ezgl::line_cap::round);  // Round ends
            g->set_line_dash(ezgl::line_dash::none); // Solid line
            g->draw_text({950, 920 + offsetY}, "Round ends, 33% transparent", 400, DBL_MAX);
        }

        else
        {
            g->set_color(ezgl::RED, 255 / 3);                  // Red line that is 67% transparent
            g->set_line_cap(ezgl::line_cap::butt);             // butt ends
            g->set_line_dash(ezgl::line_dash::asymmetric_5_3); // Dashed line
            g->draw_text({950, 920 + offsetY}, "Butt ends, 67% transparent", 400, DBL_MAX);
        }

        g->set_horiz_justification(ezgl::justification::center);

        g->draw_text({200, 900 + offsetY}, "Thin line (width 1)", 200, DBL_MAX);
        g->set_line_width(1);
        g->draw_line({100, 920 + offsetY}, {300, 920 + offsetY});

        g->draw_text({500, 900 + offsetY}, "Width 3 Line", 200, DBL_MAX);
        g->set_line_width(3);
        g->draw_line({400, 920 + offsetY}, {600, 920 + offsetY});

        g->draw_text({800, 900 + offsetY}, "Width 6 Line", 200, DBL_MAX);
        g->set_line_width(6);
        g->draw_line({700, 920 + offsetY}, {900, 920 + offsetY});

        g->set_line_width(1);
    }
}

/**
 * Draw to screen coordinates where (0,0) is the top-left corner of the window
 * These coordinates are not transformed so the object will not pan or zoom.
 */
void screen_coordinates_example(ezgl::renderer *g)
{
    // Set the coordinate system to SCREEN
    g->set_coordinate_system(ezgl::SCREEN);

    g->set_color(255, 0, 0, 255);
    g->set_line_dash(ezgl::line_dash::none);
    g->draw_rectangle({10, 10}, {100, 100});
    g->set_font_size(10);
    g->draw_text({55, 33}, "Screen coord");
    g->draw_text({55, 66}, "Fixed loc");

    // Set the coordinate system back to WORLD
    g->set_coordinate_system(ezgl::WORLD);
}

/**
 * Draw a small PNG
 */
void draw_png_example(ezgl::renderer *g)
{
    ezgl::surface *png_surface = ezgl::renderer::load_png("small_image.png");
    g->draw_surface(png_surface, {50, 200});
    ezgl::renderer::free_surface(png_surface);
    g->set_font_size(10);
    g->set_color(ezgl::BLACK);
    g->draw_text({50, 225}, "draw_surface", 200, DBL_MAX);
}

void draw_grid(ezgl::renderer *g)
{
    g->set_font_size(10);
    g->set_line_width(2);
    for (int i = 0; i <= 3; i++)
    {
        g->set_color(ezgl::DARK_KHAKI);
        g->draw_line({0 + 1000 * i, 0}, {0 + 1000 * i, 3000});
        g->draw_text({0 + 1000 * i, -40}, to_string(i * 10), 200, DBL_MAX);
        g->draw_line({0, 0 + 1000 * i}, {3000, 0 + 1000 * i});
        if (i > 0)
        {
            g->draw_text({-50, 0 + 1000 * i}, to_string(i * 10), 200, DBL_MAX);
        }
    }
}
void draw_blocks(ezgl::renderer *g)
{
    for (auto block : ap.blockCoordinates)
    {
        if (ap.blocktypes[block.first] == 0)
        {
            g->set_color(ezgl::CYAN);
            const ezgl::point2d start_point = {(block.second.first * 100) - 50, (block.second.second * 100) - 50};
            ezgl::rectangle blockToDraw = {start_point, 100, 100};
            g->fill_rectangle(blockToDraw);
            g->set_color(ezgl::BLACK);
            g->set_line_width(2);
            ezgl::rectangle border = {start_point, 100, 100};
            g->set_line_width(1);
            g->draw_text({(block.second.first * 100), (block.second.second * 100)}, to_string(block.first), 200, DBL_MAX);
        }
        else
        {
            g->set_color(ezgl::RED);
            const ezgl::point2d start_point = {(block.second.first * 100) - 100, (block.second.second * 100) - 50};
            ezgl::rectangle blockToDraw = {start_point, 200, 100};
            g->fill_rectangle(blockToDraw);
            g->set_color(ezgl::BLACK);
            g->set_line_width(2);
            ezgl::rectangle border = {start_point, 100, 100};
            g->set_line_width(1);
            g->draw_text({(block.second.first * 100), (block.second.second * 100)}, to_string(block.first), 200, DBL_MAX);
        }
    }
}
void draw_connections(ezgl::renderer *g)
{
    for (auto net : ap.nets)
    {
        for (int i = 0; i < net.second.size(); i++)
        {
            for (int j = 0; j < net.second.size(); j++)
            {
                pair<double, double> block1 = ap.blockCoordinates[net.second[i]];
                pair<double, double> block2 = ap.blockCoordinates[net.second[j]];
                g->set_color(ezgl::GREY_55);
                g->draw_line({block1.first * 100, block1.second * 100}, {block2.first * 100, block2.second * 100});
            }
        }
    }
}

void draw_HPWL(ezgl::renderer *g)
{
    g->set_color(ezgl::BLACK);
    g->draw_text({1500, -100}, "HPWL = " + to_string(ap.HPWL), 2000, DBL_MAX);
    if (CellDisplacement != 0)
    {
        g->draw_text({1500, -200}, "Total Cell Displacement = " + to_string(CellDisplacement), 2000, DBL_MAX);
    }
}
void draw_anchors(ezgl::renderer *g)
{
    for (auto i : ap.blockCoordinates)
    {
        pair<double, double> block1 = {i.second.first, i.second.second};
        pair<double, double> block2 = ap.blockCoordinates_new[i.first];
        g->set_color(ezgl::BLACK);
        g->draw_line({block1.first * 100, block1.second * 100}, {block2.first * 100, block2.second * 100});
    }
}
/**
 * A callback function to the Animate button. Creates an Animation in the main wundow
 */
void animate_button_cbk(GtkWidget * /*widget*/, ezgl::application *application)
{
    // Get a renderer that can be used to draw on top of the main canvas
    ezgl::renderer *g = application->get_renderer();

    // Set line attributes
    g->set_line_width(2);
    g->set_color(ezgl::RED);
    g->set_line_dash(ezgl::line_dash::asymmetric_5_3);

    ezgl::point2d start_point = {100, 0};

    // Do some animation
    for (int i = 0; i < 50; i++)
    {
        // Draw a line
        g->draw_line(start_point, start_point + ezgl::point2d(500, 10));
        start_point += {10, 20};

        // Flush the drawings done by the renderer to the screen
        application->flush_drawing();

        // Add some delay
        std::chrono::milliseconds duration(50);
        std::this_thread::sleep_for(duration);
    }
}

/**
 * A callback function to test the Test button. Changes application message when button is pressed
 */
void test_button_cbk(GtkWidget * /*widget*/, ezgl::application *application)
{
    // Update the status bar message
    application->update_message("Test Button Pressed");

    // Redraw the main canvas
    application->refresh_drawing();
}

/**
 * Callback function for the example combo box. Sets message to currently active option.
 * Function trigerred when currently selected option changes.
 */
void combo_box_cbk(GtkComboBoxText *self, ezgl::application *app)
{
    // Getting text content of combo box. This call makes a copy that we must free
    auto text = gtk_combo_box_text_get_active_text(self);
    if (!text)
    { // Returning if the combo box is currently empty (Always check to avoid errors)
        return;
    }
    else
    { // Updating message to reflect new combo box value.
        app->update_message(text);
        g_free(text); // gtk made a copy that we own; need to free.
    }
}

/**
 * Callback function for the delete combo box button. Deletes combo box.
 */
void delete_combo_box_cbk(GtkWidget *widget, ezgl::application *app)
{
    // Destroying widget. If function fails (could not find widget), changing message to reflect failure.
    if (app->destroy_widget("TestComboBox"))
    {
        app->update_message("Successfully deleted");
    }
    else
    {
        app->update_message("Already deleted or not found");
    }
}

/**
 * Callback function for the create dialog button. Creates a dialog window and connects it to the dialog_cbk function
 */
void create_dialog_button_cbk(GtkWidget * /*widget*/, ezgl::application *application)
{
    application->create_dialog_window(dialog_cbk, "Title", "THIS IS SOME TEXT. HELLO!");
}

/**
 * Callback function for the create message button. Creates a popup message
 */
void create_mssg_button_cbk(GtkWidget * /*widget*/, ezgl::application *app)
{
    app->create_popup_message("My Message", "Hello, hit Done to Proceed");
}

/**
 * Callback function for dialog window created by "Create Dialog Window" button.
 * Updates application message to reflect user answer to dialog window.
 */
void dialog_cbk(GtkDialog *self, gint response_id, ezgl::application *app)
{
    // Response_id is an integer/enumeration, so we can use a switch to read its value and act accordingly
    switch (response_id)
    {
    case GTK_RESPONSE_ACCEPT:
        app->update_message("USER ACCEPTED");
        break;
    case GTK_RESPONSE_REJECT:
        app->update_message("USER REJECTED");
        break;
    case GTK_RESPONSE_DELETE_EVENT:
        app->update_message("USER CLOSED WINDOW");
        break;
    default:
        app->update_message("YOU SHOULD NOT SEE THIS");
    }

    // We always have to destroy the dialog window in the callback function or it will never close
    gtk_widget_destroy(GTK_WIDGET(self));
}

/**
 * Function to handle mouse press event
 * The current mouse position in the main canvas' world coordinate system is returned
 * A pointer to the application and the entire GDK event are also returned
 */
void act_on_mouse_press(ezgl::application *application, GdkEventButton *event, double x, double y)
{
    application->update_message("Mouse Clicked");

    std::cout << "User clicked the ";

    if (event->button == 1)
        std::cout << "left ";
    else if (event->button == 2)
        std::cout << "middle ";
    else if (event->button == 3)
        std::cout << "right ";

    std::cout << "mouse button at coordinates (" << x << "," << y << ") ";

    if ((event->state & GDK_CONTROL_MASK) && (event->state & GDK_SHIFT_MASK))
        std::cout << "with control and shift pressed ";
    else if (event->state & GDK_CONTROL_MASK)
        std::cout << "with control pressed ";
    else if (event->state & GDK_SHIFT_MASK)
        std::cout << "with shift pressed ";

    std::cout
        << std::endl;
}

/**
 * Function to handle mouse move event
 * The current mouse position in the main canvas' world coordinate system is returned
 * A pointer to the application and the entire GDK event are also returned
 */
void act_on_mouse_move(ezgl::application * /*application*/, GdkEventButton * /*event*/, double x, double y)
{
    std::cout << "Mouse move at coordinates (" << x << "," << y << ") " << std::endl;
}

/**
 * Function to handle keyboard press event
 * The name of the key pressed is returned (0-9, a-z, A-Z, Up, Down, Left, Right, Shift_R, Control_L, space, Tab, ...)
 * A pointer to the application and the entire GDK event are also returned
 */
void act_on_key_press(ezgl::application *application, GdkEventKey * /*event*/, char *key_name)
{
    application->update_message("Key Pressed");

    std::cout << key_name << " key is pressed" << std::endl;
}
