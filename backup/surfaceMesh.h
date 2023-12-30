#include "IO.h"
using Eigen::Vector3d;

class SurfaceMesh{
	std::vector<Vector3d> positions;
	std::vector<Vector3d> normals;
	std::vector<uint> indices;

	std::vector<Vector3d> faceNormals;

	void computeNormals(std::vector<Vector3d> &normals_) const
	{
		normals_.resize(positions.size());
		std::fill(normals_.begin(), normals_.end(), Vector3d::Zero());

		const bool dirty = faceNormals.size() * 3 != indices.size();

		for (size_t i = 0; i < indices.size(); i += 3) {
			const Vector3d n = dirty ? getFaceNormal(i) : faceNormals[i / 3];
			for (size_t j = i; j < i + 3; j++) normals_[indices[j]] += n;
		}
		for (auto &normal : normals_) normal.normalize();
	}

	Vector3d SurfaceMesh::getFaceNormal(const size_t i) const
	{
		const Vector3d a = positions[indices[i + 1]] - positions[indices[i]];
		const Vector3d b = positions[indices[i + 2]] - positions[indices[i]];
		return a.cross(b).normalized();
	}

public:
	void write(std::ostream &out) const{
		IO::writeValue(out, uint(positions.size()));
		for (const auto &pos : positions)
		IO::writeValue(out, pos.template cast<float>().eval());
		if (normals.size() != positions.size()) {
			std::vector<Vector3d> normals_;
			computeNormals(normals_);
			for (const auto &normal : normals_)
				IO::writeValue(out, normal.template cast<float>().eval());
		}
		else {
			for (const auto &normal : normals)
				IO::writeValue(out, normal.template cast<float>().eval());
		}
		IO::writeValue(out, uint(indices.size()));
		IO::writeArray(out, indices.data(), indices.size());
	}
	
};
