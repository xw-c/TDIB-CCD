#include "Simulator.h"
#include "Simulation.h"
#include "recBezier.h"
#include "IO.h"
#include <string>
class BezierVisualize: public Simulation{
	int cntPatches;
	std::vector<RecCubicBezier> patches;
public:
	BezierVisualize(const std::string& filename){
		std::ifstream readin(filename);
		readin>>cntPatches;
		patches.resize(cntPatches);
		int uOrder, vOrder;//暂时没管
		for(auto& patch: patches){
			readin>>uOrder>>vOrder;
			for(auto& pt:patch.ctrlp)
				readin>>pt[0]>>pt[1]>>pt[2];
		}
		readin.close();
	};

	virtual ~BezierVisualize() = default;

	virtual int dimension() const override { return 3; }
	virtual void writeDescription(YAML::Node &root) const override;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string &frameDir) const override;
	virtual void loadFrame(const std::string &frameDir) override;

	// virtual void initialize() override {};
	virtual void advance(const double dt) override {};

};
void BezierVisualize::writeDescription(YAML::Node &root) const
{
	{ // Description of front.
		YAML::Node node;
		node["name"] = "front";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "triangle_list";
		node["indexed"] = true;
		node["material"]["diffuse_albedo"] = Vector4f(1, 1, 0, 1);
		root["objects"].push_back(node);
	}
	{ // Description of back.
		YAML::Node node;
		node["name"] = "back";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "triangle_list";
		node["indexed"] = true;
		node["material"]["diffuse_albedo"] = Vector4f(0, 1, 0, 1);
		root["objects"].push_back(node);
	}
	{ // Description of nodes.
		YAML::Node node;
		node["name"] = "nodes";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["indexed"] = false;
		node["material"]["diffuse_albedo"] = Vector4f(1, 0, 1, 1);
		root["objects"].push_back(node);
	}
}

void BezierVisualize::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	{ // Write patches.
		std::ofstream fout(frameDir + "/patches.dat", std::ios::binary);
		writeBezierPatches(fout);
	}
	{ // Write front and back.
		SurfaceMesh mesh = GridMesher(_grid).mesh();
		{ // Write front.
			std::ofstream fout(frameDir + "/front.mesh", std::ios::binary);
			mesh.write(fout);
		}
		{ // Write back.
			std::ofstream fout(frameDir + "/back.mesh", std::ios::binary);
			std::reverse(mesh.indices.begin(), mesh.indices.end());
			mesh.write(fout);
		}
	}
	{ // Write nodes.
		std::ofstream fout(frameDir + "/nodes.mesh", std::ios::binary);
		IO::writeValue(fout, uint(_grid.nodeCnt));
		_grid.forEachNode([&](const int i, const int j) {
			IO::writeValue(fout, _grid(i, j)[0].val.cast<float>().eval());
		});
		_grid.forEachNode([&](const int i, const int j) {
			const Vector3d normal = _grid(i, j)[1].val.cross(_grid(i, j)[2].val).normalized();
			IO::writeValue(fout, normal.cast<float>().eval());
		});
	}
}

void BezierVisualize::saveFrame(const std::string &frameDir) const
{
	{ // Save nodeDoFs.
		std::ofstream fout(frameDir + "/nodeDoFs.sav", std::ios::binary);
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp)
				IO::writeValue(fout, pt);
	}
}

void BezierVisualize::loadFrame(const std::string &frameDir)
{
	initialize();
	{ // Load nodeDoFs.
		std::ifstream fin(frameDir + "/nodeDoFs.sav", std::ios::binary);
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp)
				IO::readValue(fin, pt);
	}
}
int main(int argc, char *argv[])
{
	const std::string output ("output");
	const uint begin = 0, end = 1, rate = 50, step = 1;

	auto simulator = std::make_unique<Simulator>(output, begin, end, rate, step, BezierVisualize.get());
	simulator->Simulate();

	return 0;
}