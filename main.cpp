//212-Кончугаров-Тимур
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <cstdint>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <random>
#include <set>
#include <stack> 

//212-Кончугаров-Тимур
class Server { 
    std::ofstream server_log; 
    bool logging_enabled = false; 
    std::string log_filename = "server_log.txt"; 
public:
    Server() {
    }

    ~Server() {
        if (server_log.is_open()) {
            server_log.close();
        }
    }

    void enable_logging(const std::string& filename) { 
        logging_enabled = true;
        log_filename = filename;
        if (server_log.is_open()) {
            server_log.close();
        }
        server_log.open(log_filename, std::ios::app);
    }

    void disable_logging() {
        logging_enabled = false;
        if (server_log.is_open()) {
            server_log.close();
        }
    }

    void log(const std::string& message) { 
        if (logging_enabled && server_log.is_open()) {
            server_log << current_time() << " - " << message << std::endl;
        }
    }

    std::string current_time() { 
        std::time_t t = std::time(nullptr);
        std::tm* now = std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(now, "%d.%m.%Y %H:%M:%S");
        return oss.str();
    }
};

//212-Кончугаров-Тимур
class Client { 
    std::ofstream client_log;  
    bool logging_enabled = false; 
    std::string log_filename = "client_log.txt"; 

public:
    Client() {
    }

    ~Client() {
        if (client_log.is_open()) {
            client_log.close();
        }
    }

    void enable_logging(const std::string& filename) {
        logging_enabled = true;
        log_filename = filename;
        if (client_log.is_open()) {
            client_log.close();
        }
        client_log.open(log_filename, std::ios::app);
    }

    void disable_logging() { 
        logging_enabled = false;
        if (client_log.is_open()) {
            client_log.close();
        }
    }

    void log(const std::string& message) {
        if (logging_enabled && client_log.is_open()) {
            client_log << current_time() << " - " << message << std::endl;
        }
    }

    std::string current_time() { 
        std::time_t t = std::time(nullptr);
        std::tm* now = std::localtime(&t);
        std::ostringstream oss;
        oss << std::put_time(now, "%d.%m.%Y %H:%M:%S");
        return oss.str();
    }
};

struct Point {
    int x, y;
};

struct Color {
    uint8_t r, g, b;
};

//212-Кончугаров-Тимур
class Component {
public:
    std::vector<std::vector<double>> comp;
    Component(const std::vector<std::vector<double>>& inpcomp) : comp(inpcomp) {}

    Component(int A, int B) {
        comp.resize(A, std::vector<double>(B, 255));
    }
};


class Gauss { 
public:
    double h, x0, y0, sigma_x, sigma_y;
    Gauss(double h, double x0, double y0, double sigma_x, double sigma_y)
        : h(h), x0(x0), y0(y0), sigma_x(sigma_x), sigma_y(sigma_y) {}

    double calculate(int x, int y) const { 
        double dx = x - x0;
        double dy = y - y0;
        double denom_x = 2 * sigma_x * sigma_x;
        double denom_y = 2 * sigma_y * sigma_y;

        if ((int)denom_x == 0 || (int)denom_y == 0) {
            std::cerr << "Error: sigma_x or sigma_y iz zero!" << std::endl;
            return 0;
        }

        double exponent = -((dx * dx) / denom_x +
            (dy * dy) / denom_y);

        if (exponent < -700) {
            return 0;
        }

        else if (exponent > 700) {
            return std::numeric_limits<double>::infinity();
        }

        return h * exp(exponent);
    }
};

//212-Кончугаров-Тимур
class Field { 
public:
    int length, width;
    std::vector<std::vector<double>> matrix;
    std::vector<std::vector<double>> weight_matrix;
    std::vector<std::vector<Color>> colors; 

    Field(int l, int w) : length(l), width(w) {
        matrix.resize(length, std::vector<double>(width, 0));
        weight_matrix.resize(length, std::vector<double>(width, 0));
    }

    void apply_gauss(const Gauss& g) {
        double value;
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < width; ++j) {
                value = g.calculate(i, j);
                matrix[i][j] += value * g.h; 
                weight_matrix[i][j] += g.h; 
            }
        }
    }

    void normalize(void) {
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < width; ++j) {
                if (weight_matrix[i][j] > 0) {
                    matrix[i][j] /= weight_matrix[i][j];
                }
            }
        }
    }

    void save_to_gnuplot(const std::string& filename) { 
        std::ofstream bmp_file(filename);
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < width; ++j) {
                bmp_file << i << " " << j << " " << matrix[i][j] << std::endl;
            }
            bmp_file << std::endl;
        }
        bmp_file.close();
    }

    void bmp_write(const std::string& filename, int k) { 
        const int BMP_HEADER_SIZE = 54;
        const int PIXEL_SIZE = 3;
        int file_size = BMP_HEADER_SIZE + PIXEL_SIZE * length * width;

        unsigned char bmp_header[BMP_HEADER_SIZE] = { 0 };

        bmp_header[0] = 'B';
        bmp_header[1] = 'M';
        bmp_header[2] = file_size & 0xFF;
        bmp_header[3] = (file_size >> 8) & 0xFF;
        bmp_header[4] = (file_size >> 16) & 0xFF;
        bmp_header[5] = (file_size >> 24) & 0xFF;
        bmp_header[10] = BMP_HEADER_SIZE;

        bmp_header[14] = 40; 
        bmp_header[18] = width & 0xFF;
        bmp_header[19] = (width >> 8) & 0xFF;
        bmp_header[20] = (width >> 16) & 0xFF;
        bmp_header[21] = (width >> 24) & 0xFF;
        bmp_header[22] = length & 0xFF;
        bmp_header[23] = (length >> 8) & 0xFF;
        bmp_header[24] = (length >> 16) & 0xFF;
        bmp_header[25] = (length >> 24) & 0xFF;
        bmp_header[26] = 1; 
        bmp_header[28] = 24; 

        std::ofstream bmp_file(filename, std::ios::binary);
        bmp_file.write(reinterpret_cast<char*>(bmp_header), BMP_HEADER_SIZE);

        if (k == 1) { 
            for (int i = length - 1; i >= 0; --i) {
                for (int j = 0; j < width; ++j) {
                    int value = static_cast<int>(matrix[i][j] * 100); 
                    unsigned char pixel[3] = { static_cast<unsigned char>(std::min(std::max(value, 0), 255)),
                                            static_cast<unsigned char>(std::min(std::max(value, 0), 255)),
                                            static_cast<unsigned char>(std::min(std::max(value, 0), 255)) };
                    bmp_file.write(reinterpret_cast<char*>(pixel), PIXEL_SIZE);
                }
            }
        }

        else { 
            for (int i = length - 1; i >= 0; --i) {
                for (int j = 0; j < width; ++j) {
                    unsigned char pixel[3] = { static_cast<unsigned char>(colors[i][j].r),
                                               static_cast<unsigned char>(colors[i][j].g),
                                               static_cast<unsigned char>(colors[i][j].b) };
                    bmp_file.write(reinterpret_cast<char*>(pixel), PIXEL_SIZE);
                }
            }
        }

        bmp_file.close();
    }

    std::pair<int, int> bmp_read(const std::string& filename) { 

        std::ifstream bmp_file(filename, std::ios::binary);

        if (!bmp_file) {
            std::cerr << "Failed to open BMP file: " << filename << std::endl;
            return { 0,0 };
        }

        unsigned char header[54];
        bmp_file.read(reinterpret_cast<char*>(header), 54);

        if (header[0] != 'B' || header[1] != 'M') {
            std::cerr << "Invalid BMP file:" << filename << std::endl;
            return { 0,0 };
        }

        int width = header[18] | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
        int height = header[22] | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);

        return { height, width };
    }

    void load_data(std::ifstream& bmp_file, int length, int width) { 
        for (int i = length - 1; i >= 0; --i) {
            for (int j = 0; j < width; j++) {
                unsigned char color = bmp_file.get();
                bmp_file.get();
                bmp_file.get();
                matrix[i][j] = color;
            }
            bmp_file.ignore((4 - (width * 3) % 4) % 4);
        }
    }
};

//212-Кончугаров-Тимур
class Razrez {
    int length, width;
    std::vector<std::vector<double>> CopyField;
    std::vector<Component> components;
    int count = 0;
public:

    void bin(int h, Field& f) { 
        CopyField.resize(f.matrix.size(), std::vector<double>(f.matrix[0].size(), 0));
        length = f.matrix.size();
        width = f.matrix[0].size();

        for (int i = length - 1; i >= 0; --i) {
            for (int j = 0; j < width; ++j) {
                if (f.matrix[i][j] < h) {
                    CopyField[i][j] = 255;
                }
                else {
                    CopyField[i][j] = 0;
                }
            }
        }
        Field pole(f.length, f.width);
        pole.matrix = CopyField;
        pole.bmp_write("bin.bmp", 1);
    }

    void wave(int n) {
        Component Componenta(length, width);
        int c = 0; 

        for (int i = length - 1; i >= 0; --i) {
            for (int j = 0; j < width; ++j) {
                c = inc(Componenta.comp, i, j, 1);

                if (c > n) {
                    components.emplace_back(Componenta);
                    Componenta = Component(length, width);
                }

                else if (c > 0 && c < n) {
                    Componenta = Component(length, width);
                }

                count = 0;

            }
        }

        for (int i = 0; i < (int)components.size(); i++) {
            Field compole(length, width);
            compole.matrix = components[i].comp;
            compole.bmp_write("Comp" + std::to_string(i + 1) + ".bmp", 1);
        }
    }

    int inc(std::vector<std::vector<double>>& component, int x, int y, int k) {
        if (x < 0 || y < 0 || x >= (int)CopyField.size() || y >= (int)CopyField[0].size() || CopyField[x][y] == 255)
            return 0;

        std::stack<std::pair<int, int>> stack; 
        stack.push({ x, y });
        int count = 0;

        while (!stack.empty()) {
            std::pair<int, int> current = stack.top();
            stack.pop();
            int cx = current.first;
            int cy = current.second;

            if (cx < 0 || cy < 0 || cx >= (int)CopyField.size() || cy >= (int)CopyField[0].size() || CopyField[cx][cy] == 255)
                continue;

            CopyField[cx][cy] = 255;
            component[cx][cy] = 0;
            count++;

            stack.push({ cx + 1, cy });
            stack.push({ cx - 1, cy });
            stack.push({ cx, cy + 1 });
            stack.push({ cx, cy - 1 });
            stack.push({ cx + 1, cy + 1 });
            stack.push({ cx - 1, cy - 1 });
            stack.push({ cx + 1, cy - 1 });
            stack.push({ cx - 1, cy + 1 });
        }

        return count;
    }

    std::vector<Color> generateColor(int k) {
        std::vector<Color> colors;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(0, 255);

        for (int i = 0; i < k; i++) {
            colors.push_back({ static_cast<uint8_t>(dis(gen)),
                              static_cast<uint8_t>(dis(gen)),
                              static_cast<uint8_t>(dis(gen)) });
        }

        return colors;
    }

    void kMeans(int k, int p) {

        std::vector<Point> points;

        for (const auto& comp : components) {
            for (int i = 0; i < length; i++) {
                for (int j = 0; j < width; j++) {
                    if (comp.comp[i][j] == 0) {
                        points.push_back({ i, j });
                    }
                }
            }
        }

        if (points.empty()) {
            std::cerr << "No points available for clustering!" << std::endl;
            return;
        }

        std::vector<Point> centroids;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis_x(0, length - 1);
        std::uniform_int_distribution<> dis_y(0, width - 1);

        for (int i = 0; i < k; i++) {
            centroids.push_back({ dis_x(gen), dis_y(gen) }); // random pic k centroids
        }

        bool changed = true;
        std::vector<int> labels(points.size(), -1);

        while (changed) {
            changed = false;

            //1: pic point for nearest centroids
            for (size_t i = 0; i < points.size(); i++) {
                double minDist = std::numeric_limits<double>::max();
                int cluster = -1;

                for (int j = 0; j < k; j++) { //count lenght from point to cluster, pic min
                    double dist = std::pow(points[i].x - centroids[j].x, 2) + std::pow(points[i].y - centroids[j].y, 2);
                    if (dist < minDist) {
                        minDist = dist;
                        cluster = j;
                    }
                }

                if (labels[i] != cluster) { //check if point change cluster
                    labels[i] = cluster; // if change, change her cluster
                    changed = true;
                }
            }

            //2: calculate new centroids
            std::vector<Point> newCentroids(k, { 0, 0 });
            std::vector<int> counts(k, 0);

            for (size_t i = 0; i < points.size(); i++) { 
                newCentroids[labels[i]].x += points[i].x;
                newCentroids[labels[i]].y += points[i].y;
                counts[labels[i]]++;
            }

            for (int j = 0; j < k; j++) {
                if (counts[j] > 0) {
                    newCentroids[j].x /= counts[j];
                    newCentroids[j].y /= counts[j];
                }
            }

            centroids = newCentroids;

        }

        std::vector<Color> clusterColors = generateColor(k); //generate colors for clusters
        std::vector<std::vector<Color>> clusterImage(length, std::vector<Color>(width, { 255, 255, 255 })); //Заополняем матрицу цветов каждого пикселя
        std::vector<std::vector<Point>> clusteredPoints(k);

        for (size_t i = 0; i < points.size(); i++) { //point - сщдщк
            int cluster = labels[i];
            clusterImage[points[i].x][points[i].y] = clusterColors[cluster];
            clusteredPoints[cluster].push_back(points[i]); //+every point in clusters
        }

        // 2: kMeans for every k clusters with p centroids
        std::vector<std::vector<Point>> subCentroids(k);

        for (int clusterIdx = 0; clusterIdx < k; clusterIdx++) { //kMeans for every clusters
            if (clusteredPoints[clusterIdx].empty()) continue;

            std::vector<Point>& clusterPoints = clusteredPoints[clusterIdx];
            std::vector<Point> subCenters(p);

            //init p random subcentroids
            for (int i = 0; i < p; i++) {
                subCenters[i] = clusterPoints[std::uniform_int_distribution<>(0, clusterPoints.size() - 1)(gen)];
            }

            bool subChanged = true;
            std::vector<int> subLabels(clusterPoints.size(), -1);

            while (subChanged) {
                subChanged = false;

                //point 
                for (size_t i = 0; i < clusterPoints.size(); i++) {
                    double minDist = std::numeric_limits<double>::max();
                    int subCluster = -1;

                    for (int j = 0; j < p; j++) {
                        double dist = std::pow(clusterPoints[i].x - subCenters[j].x, 2) +
                            std::pow(clusterPoints[i].y - subCenters[j].y, 2);
                        if (dist < minDist) {
                            minDist = dist;
                            subCluster = j;
                        }
                    }

                    if (subLabels[i] != subCluster) {
                        subLabels[i] = subCluster;
                        subChanged = true;
                    }
                }

                //count subcentroids
                std::vector<Point> newSubCenters(p, { 0, 0 });
                std::vector<int> subCounts(p, 0);

                for (size_t i = 0; i < clusterPoints.size(); i++) {
                    newSubCenters[subLabels[i]].x += clusterPoints[i].x;
                    newSubCenters[subLabels[i]].y += clusterPoints[i].y;
                    subCounts[subLabels[i]]++;
                }

                for (int j = 0; j < p; j++) {
                    if (subCounts[j] != 0) {
                        newSubCenters[j].x /= subCounts[j];
                        newSubCenters[j].y /= subCounts[j];
                    }
                }

                subCenters = newSubCenters;
            }
            subCentroids[clusterIdx] = subCenters;

            for (const auto& subCenter : subCenters) {//pict subcentroids
                clusterImage[subCenter.x][subCenter.y] = { 0, 0, 0 };
            }
        }

        Field pole(length, width);
        pole.colors = clusterImage;
        pole.bmp_write("kmeans.bmp", 2);
    }
};

//212-Кончугаров-Тимур
class Control {
    std::vector<Gauss> gausses;
    Field* field = nullptr;
    Razrez razrez;

public:
    Server server;

    void wave_cntrl(int n) {

        if (!field) {
            std::cerr << "Error! Field not initialized!" << std::endl;
            server.log("Control done 'bin', but field was not initialized.");
            return;
        }

        razrez.wave(n);
        server.log("Control completed the wave program");
    }

    void init(int length, int width) { 
        if (field) {
            server.log("Control done 'init', but field was already initialized.");
            return;
        }
        field = new Field(length, width);
        server.log("Field done " + std::to_string(length) + "x" + std::to_string(width));
        server.log("Control done 'init' command to Field");
    }

    void bin_cntrl(int h) {
        if (!field) {
            std::cerr << "Error! Field not initialized!" << std::endl;
            server.log("Control done 'bin', but field was not initialized.");
            return;
        }
        razrez.bin(h, *field);
        server.log("A section on height " + std::to_string(h) + " is obtained");
    }

    void kMeans_cntrl(int k, int p) {
        if (!field) {
            std::cerr << "Error! Field not initialized!" << std::endl;
            server.log("Control done 'KMEANS', but field was not initialized.");
            return;
        }
        server.log("KMEANS for the k =" + std::to_string(k) + " is launched");
        razrez.kMeans(k, p);
        server.log("KMEANS for the k =" + std::to_string(k) + " is completed");
    }

    void add_gauss(double h, double x0, double y0, double sigma_x, double sigma_y) {

        if (!field) {
            std::cerr << "Error! Field not initialized!" << std::endl;
            server.log("Control done 'g', but field was not initialized.");
            return;
        }

        gausses.emplace_back(h, x0, y0, sigma_x, sigma_y);

        server.log("Control done Gauss command to Field with (h=" + std::to_string(h) + ", x0=" + std::to_string(x0) +
            ", y0=" + std::to_string(y0) + ", sigma_x=" + std::to_string(sigma_x) + ", sigma_y=" + std::to_string(sigma_y) + ")");
    }


    void generate() {
        if (!field) {
            std::cerr << "Error! Field not initialized!" << std::endl;
            server.log("Control done 'generate', but field was not initialized.");
            return;
        }

        if (gausses.empty()) {
            server.log("No gausses to apply in 'generate'");
            return;
        }

        server.log("Control done 'generate' command to Field");

        for (size_t i = 0; i < gausses.size(); ++i) {
            server.log("Control is applying Gauss #" + std::to_string(i + 1) + " to Field");
            field->apply_gauss(gausses[i]);
            server.log("Gauss #" + std::to_string(i + 1) + "applied");
        }

        field->normalize();
        gausses.clear();
        server.log("Gauss completed applying all Gausses");
    }


    void gnuplot(const std::string& filename) {
        if (!field) {
            std::cerr << "Error! Field not initialized!" << std::endl;
            server.log("Control done 'gnuplot', but field was not initialized.");
            return;
        }

        field->save_to_gnuplot("gnuplot_" + filename + ".txt");
        server.log("Field saved to Gnuplot file: " + filename);

        std::ofstream gp_file("gnuplot_commands.txt");
        gp_file << "set view 60,30\n";
        gp_file << "set palette defined (0 \"blue\", 1 \"red\")\n";
        gp_file << "set pm3d at s\n";
        gp_file << "splot 'gnuplot_" << filename << ".txt' with lines\n";
        gp_file << "pause -1";
        gp_file.close();
        server.log("Field generated Gnuplot file: " + filename);
    }

    void bmp_write_cntrl(const std::string& filename, int k) {
        if (!field) {
            std::cerr << "Error! Field not initialized!" << std::endl;
            server.log("Control done 'save bmp', but field was not initialized.");
            return;
        }
        field->bmp_write(filename, k);
        server.log("Field saved to BMP file:" + filename);
    }

    void bmp_read_cntrl(const std::string& filename) {
        if (!field) {
            std::cerr << "Error! Field not initialized!" << std::endl;
            server.log("Control done 'read bmp', but field was not initialized.");
            return;
        }

        std::pair<int, int> dimensions = field->bmp_read(filename);
        int new_length = dimensions.first;
        int new_width = dimensions.second;

        if (new_length == 0 || new_width == 0) {
            std::cerr << "Error reading BMP file. No changes made" << std::endl;
            server.log("Error reading BMP file");
        }

        delete field;
        field = new Field(new_length, new_width);

        std::ifstream bmp_file(filename, std::ios::binary);
        bmp_file.ignore(54);
        field->load_data(bmp_file, new_length, new_width);
        bmp_file.close();

        server.log("Field loaded bmp file:" + filename);
    }
};

//212-Кончугаров-Тимур
class Interface {
    Control control;
    Client client;

public:

    std::string trim(const std::string& str) {
        size_t first = str.find_first_not_of(' ');
        if (first == std::string::npos) {
            return "";
        }
        size_t last = str.find_last_not_of(' ');
        return str.substr(first, (last - first + 1));
    }

    int run() {
        std::string filename;
        std::cout << "Let's open the config file: ";
        std::cin >> filename;

        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  
        std::ifstream config(filename);
        bool gaussi = false;
        bool batch_pause = false;
        int n = 0; 

        if (!config.is_open()) {
            std::cerr << "Failed to open file:" << filename << std::endl;
            return 1;
        }

        std::cout << "Config file opened." << std::endl;
        std::string line;
        bool noize = false;

        while (!noize && std::getline(config, line)) {
            std::istringstream iss(line);
            std::string key, value;

            if (std::getline(iss, key, '=') && std::getline(iss, value)) {
                key = trim(key);
                value = trim(value);

                if (key == "Noize") {
                    std::istringstream ss(value);
                    ss >> n;
                    noize = true;
                }
            }
        }

        config.clear();
        config.seekg(0, std::ios::beg);

        while (std::getline(config, line)) {
            std::istringstream iss(line);
            std::string key, value;

            if (std::getline(iss, key, '=') && std::getline(iss, value)) {
                key = trim(key);
                value = trim(value);

                if (key == "Log_Server") {
                    if (value == "ON") {
                        std::string log_filename = "server_log.txt";
                        std::getline(config, line);
                        std::istringstream iss_filename(line);
                        if (std::getline(iss_filename, key, '=') && std::getline(iss_filename, value)) {
                            key = trim(key);
                            value = trim(value);
                            if (key == "Log_Server_filename") {
                                log_filename = value;
                            }
                        }
                        control.server.enable_logging(log_filename);
                        std::cout << "Server logs are written to the file: " << log_filename << std::endl;
                    }
                }


                else if (key == "Log_Client") {

                    if (value == "ON") {
                        std::string log_filename = "client_log.txt";
                        std::getline(config, line);
                        std::istringstream iss_filename(line);

                        if (std::getline(iss_filename, key, '=') && std::getline(iss_filename, value)) {
                            key = trim(key);
                            value = trim(value);
                            if (key == "Log_Client_filename") {
                                log_filename = value;
                            }
                        }

                        client.enable_logging(log_filename);
                        std::cout << "Client logs are written to the file: " << log_filename << std::endl;
                        /*
                        if (!noize) {
                            client.log("The noise value was not found, so we assumed that there is no noise");
                        }

                        else {
                            client.log("The noise value is set to" + std::to_string(n));
                        }*/
                    }
                }

                else if (key == "Batch_pause") {
                    if (value == "ON") {
                        batch_pause = true;
                    }
                    else {
                        batch_pause = false;
                    }
                }

                else if (key == "Use_file") {
                    if (value == "ON") {

                        std::string Commands_file = "Use_file.txt";
                        std::getline(config, line);
                        std::istringstream iss_filename(line);

                        if (std::getline(iss_filename, key, '=') && std::getline(iss_filename, value)) {
                            key = trim(key);
                            value = trim(value);
                            if (key == "Commands_file") {
                                Commands_file = value;
                            }
                        }

                        std::ifstream batch(Commands_file);

                        if (!batch.is_open()) {
                            std::cerr << "Failed to open file:" << Commands_file << std::endl;
                            return 1;
                        }

                        std::cout << "The commands are taken from the file: " << Commands_file << std::endl;

                        std::string batch_line, batch_key, batch_value;
                        std::getline(batch, batch_line);
                        std::istringstream iss_batch(batch_line);
                        iss_batch >> batch_key;

                        if (batch_key != "init") {
                            std::getline(config, line);
                            iss.clear();
                            iss.str(line);
                            std::getline(iss, key, '=');
                            std::getline(iss, value);
                            key = trim(key);
                            value = trim(value);

                            if (key == "Default_field") {
                                int length, width;
                                std::istringstream iss_Field(value);
                                iss_Field >> length >> width;
                                client.log("The interface has set the default value of the field: length=" + std::to_string(length) + ", width=" + std::to_string(width));
                                control.init(length, width);
                            }

                            else {
                                client.log("Field initialization error");
                                return -1;
                            }

                            batch.clear(); 
                            batch.seekg(0, std::ios::beg);
                        }

                        else {
                            batch.clear();
                            batch.seekg(0, std::ios::beg);
                        }


                        while (std::getline(batch, batch_line)) {

                            if (batch_pause == true) {
                                std::cout << "Нажмите Enter, чтобы продолжить.." << std::endl;
                                std::cin.get();
                            }
                            std::istringstream iss_batch(batch_line);

                            if (iss_batch >> batch_key) {
                                std::getline(iss_batch, batch_value);
                                std::istringstream iss_info(batch_value);

                                if (batch_key == "init") {
                                    int length, width;
                                    iss_info >> length >> width;
                                    if (batch_pause == true) {
                                        std::cout << "Interface done 'init' command: length=" + std::to_string(length) + ", width=" + std::to_string(width) << std::endl;
                                    }
                                    client.log("Interface done 'init' command: length=" + std::to_string(length) + ", width=" + std::to_string(width));
                                    control.init(length, width);

                                }
                                else if (batch_key == "g") {
                                    double h, x0, y0, sigma_x, sigma_y;
                                    iss_info >> h >> x0 >> y0 >> sigma_x >> sigma_y;
                                    if (batch_pause == true) {
                                        std::cout << "Interface done 'g' command: h=" + std::to_string(h) + ", x0=" + std::to_string(x0) + ", y0=" + std::to_string(y0) +
                                            ", sigma_x=" + std::to_string(sigma_x) + ", sigma_y=" + std::to_string(sigma_y) << std::endl;
                                    }
                                    client.log("Interface done 'g' command: h=" + std::to_string(h) + ", x0=" + std::to_string(x0) + ", y0=" + std::to_string(y0) +
                                        ", sigma_x=" + std::to_string(sigma_x) + ", sigma_y=" + std::to_string(sigma_y));
                                    control.add_gauss(h, x0, y0, sigma_x, sigma_y);
                                    gaussi = true;

                                }
                                else if (batch_key == "generate") {
                                    if (gaussi == false) {
                                        std::streampos saved_position = config.tellg();
                                        while (std::getline(config, line)) {
                                            iss.str(line);
                                            iss.clear();
                                            if (std::getline(iss, key, '=') && std::getline(iss, value)) {
                                                key = trim(key);
                                                value = trim(value);
                                                if (key == "Default_gauss") {
                                                    double h, x0, y0, sigma_x, sigma_y;
                                                    iss.clear();
                                                    iss.str(value);
                                                    iss >> h >> x0 >> y0 >> sigma_x >> sigma_y;
                                                    if (batch_pause == true) {
                                                        std::cout << "The default Gauss value is set: h=" + std::to_string(h) + ", x0=" + std::to_string(x0) + ", y0=" + std::to_string(y0) +
                                                            ", sigma_x=" + std::to_string(sigma_x) + ", sigma_y=" + std::to_string(sigma_y) << std::endl;
                                                    }
                                                    client.log("The default Gauss value is set: h=" + std::to_string(h) + ", x0=" + std::to_string(x0) + ", y0=" + std::to_string(y0) +
                                                        ", sigma_x=" + std::to_string(sigma_x) + ", sigma_y=" + std::to_string(sigma_y));
                                                    control.add_gauss(h, x0, y0, sigma_x, sigma_y);
                                                    gaussi = true;
                                                }
                                            }
                                        }
                                        config.seekg(saved_position);
                                    }

                                    if (batch_pause == true) {
                                        std::cout << "Interface done 'generate' command." << std::endl;
                                    }
                                    client.log("Interface done 'generate' command.");
                                    control.generate();

                                }
                                else if (batch_key == "gnuplot") {
                                    std::string filename;
                                    iss_info >> filename;

                                    if (batch_pause == true) {
                                        std::cout << "Interface done 'gnuplot' command with filename: " + filename << std::endl;
                                    }

                                    client.log("Interface done 'gnuplot' command with filename: " + filename);
                                    control.gnuplot(filename);

                                }
                                else if (batch_key == "bmp") {
                                    std::string filename;
                                    iss_info >> filename;
                                    if (batch_pause == true) {
                                        std::cout << "Interface done 'write bmp' command with filename:" + filename << std::endl;
                                    }
                                    client.log("Interface done 'write bmp' command with filename:" + filename);
                                    control.bmp_write_cntrl(filename, 1);

                                }
                                else if (batch_key == "read_bmp") {
                                    std::string filename;
                                    iss_info >> filename;
                                    if (batch_pause == true) {
                                        std::cout << "Interface done 'read bmp' command with filename: " + filename << std::endl;
                                    }
                                    client.log("Interface done 'read bmp' command with filename: " + filename);
                                    control.bmp_read_cntrl(filename);
                                    gaussi = true;

                                }
                                else if (batch_key == "bin") {
                                    int h;
                                    iss_info >> h;

                                    if (batch_pause == true) {
                                        std::cout << "Interface done 'bin' command with h = " + std::to_string(h) << std::endl;
                                    }

                                    client.log("Interface done 'bin' command with h = " + std::to_string(h));
                                    control.bin_cntrl(h);

                                }
                                else if (batch_key == "wave") {

                                    if (batch_pause == true) {
                                        std::cout << "Interface done 'wave' command with pixel noise:" + std::to_string(n) << std::endl;
                                    }

                                    client.log("Interface done 'wave' command with pixel noise:" + std::to_string(n));
                                    control.wave_cntrl(n);
                                }
                                else if (batch_key == "kmeans") {
                                    int k, p;
                                    iss_info >> k >> p;

                                    if (batch_pause == true) {
                                        std::cout << "Interface done 'kmeans' coommand with k =" + std::to_string(k) << std::endl;
                                    }

                                    client.log("Interface done 'kmeans' coommand with k =" + std::to_string(k));
                                    control.kMeans_cntrl(k, p);
                                }
                            }

                        }
                    }
                }
            }
        }
        return 0;
    }
};

//212-Кончугаров-Тимур
int main() {
    Interface interface;
    interface.run();
    return 0;
}
