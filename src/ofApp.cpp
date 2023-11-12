#include "ofApp.h"

#include "glm/gtx/euler_angles.hpp"
#include <math.h>

int counter = 0;
// Mesh functions
bool Mesh::intersect(const Ray& ray, glm::vec3& point, glm::vec3& normal) {
	glm::vec2 baricenter;
	float distance = 0;
	float shortest = std::numeric_limits<float>::max();

	bool hit = false;
	int intersectedTriangleIndex = -1;
	for (int i = 0; i < triList.size(); i++) {
		Triangle t = triList[i];
		if (glm::intersectRayTriangle(ray.p, ray.d, mVertices[t.vertIndices[0]].pos, mVertices[t.vertIndices[1]].pos, mVertices[t.vertIndices[2]].pos, baricenter, distance)) {
			if (distance > 0.0 && distance < shortest) {
				shortest = distance;
				intersectedTriangleIndex = i;
				hit = true;
			}

		}
	}

	if (hit) {
		Ray r = ray;
		point = r.evalPoint(shortest);
		normal = getTriangleNormal(intersectedTriangleIndex);
	}

	return hit;
}

void Mesh::draw() {
	//iterates through list of triangles and draw each one
	for (int i = 0; i < triList.size(); i++) {
		Triangle tri = triList[i];//current tri
		if (tri.vertIndices.size() != 3) //check if this triangle has 3 vertices
			throw std::invalid_argument("Triangle does not have exactly 3 vertices");
		else {
			glm::vec3 p1 = mVertices[tri.vertIndices[0]].pos;//position of the vertex of the index stored by Triangle
			glm::vec3 p2 = mVertices[tri.vertIndices[1]].pos;
			glm::vec3 p3 = mVertices[tri.vertIndices[2]].pos;
			ofSetColor(ofColor::white);
			ofDrawLine(p1, p2);
			ofDrawLine(p2, p3);
			ofDrawLine(p1, p3);
			//			glm::vec3 centroid = (p1 + p2 + p3) / 3;
			//			glm::vec3 normal = getTriangleNormal(i);
			//			ofDrawLine(centroid,centroid + normal);
		}
	}
}

// Intersect Ray with Plane  (wrapper on glm::intersect*
//
bool Plane::intersect(const Ray& ray, glm::vec3& point, glm::vec3& normalAtIntersect) {
	float dist;
	bool insidePlane = false;
	bool hit = glm::intersectRayPlane(ray.p, ray.d, position, this->normal, dist);
	if (hit) {
		Ray r = ray;
		point = r.evalPoint(dist);
		normalAtIntersect = this->normal;
		glm::vec2 xrange = glm::vec2(position.x - width / 2, position.x + width / 2);
		glm::vec2 yrange = glm::vec2(position.y - width / 2, position.y + width / 2);
		glm::vec2 zrange = glm::vec2(position.z - height / 2, position.z + height / 2);

		// horizontal 
		//
		if (normal == glm::vec3(0, 1, 0) || normal == glm::vec3(0, -1, 0)) {
			if (point.x < xrange[1] && point.x > xrange[0] && point.z < zrange[1] && point.z > zrange[0]) {
				insidePlane = true;
			}
		}
		// front or back
		//
		else if (normal == glm::vec3(0, 0, 1) || normal == glm::vec3(0, 0, -1)) {
			if (point.x < xrange[1] && point.x > xrange[0] && point.y < yrange[1] && point.y > yrange[0]) {
				insidePlane = true;
			}
		}
		// left or right
		//
		else if (normal == glm::vec3(1, 0, 0) || normal == glm::vec3(-1, 0, 0)) {
			if (point.y < yrange[1] && point.y > yrange[0] && point.z < zrange[1] && point.z > zrange[0]) {
				insidePlane = true;
			}
		}
	}
	return insidePlane;
}

//determine which pixel on image to map to u,v on Plane
//with texturingMode true for texture, false for specular
ofColor Plane::texture_lookup_wrap(glm::vec3 shadingPoint, bool useSpec) {
	if (!textureImage.isAllocated())
		return diffuseColor;
	ofImage* image = &textureImage;
	if (useSpec)
		image = &specularImage;
	// (x',y',z') -> (x,y) convert world to object's space
	// subtract shadingpoint by plane's position to get position relative to plane
	glm::vec3 toUse = shadingPoint - position;
	//(x,y) -> (u,v) 
	float u = ofMap(toUse.x, 0, image->getWidth(), 0, 10);
	float v = ofMap(toUse.y, 0, image->getHeight(), 0, 10);
	float w = ofMap(toUse.z, 0, image->getHeight(), 0, 10);
	//(u,v) -> (i,j)
	float temp1 = u * image->getWidth() - 0.5;
	float temp2 = 0;
	float i = temp1 - image->getWidth() * floor(temp1 / image->getWidth());
	//float i = glm::max(0.0f, glm::min(image->getWidth() - 1, temp1));
	float j = 0.0;
	if (normal.y == 1 || normal.y == -1) {//plane is floor or ceiling
		temp2 = w * image->getHeight() - 0.5;

	}
	else if (normal.z == 1) {//plane is wall
		temp2 = v * image->getHeight() - 0.5;
	}
	j = temp2 - image->getHeight() * floor(temp2 / image->getHeight());
	//j = glm::max(0.0f, glm::min(image->getHeight() - 1, temp2));
	return image->getColor(i, j);
}



// Convert (u, v) to (x, y, z) 
// We assume u,v is in [0, 1]
//
glm::vec3 ViewPlane::toWorld(float u, float v) {
	float w = width();
	float h = height();
	return (glm::vec3((u * w) + min.x, (v * h) + min.y, position.z));
}

// Get a ray from the current camera position to the (u, v) position on
// the ViewPlane
//
Ray RenderCam::getRay(float u, float v) {
	glm::vec3 pointOnPlane = view.toWorld(u, v);
	return(Ray(position, glm::normalize(pointOnPlane - position)));
}

void SpotLight::draw() {
	//draw aim direction
	ofSetColor(ofColor::green);
	ofDrawLine(position, aimPoint);

	ofSetColor(ofColor::red);
	glm::vec3 spherePos = position + (glm::normalize(aimPoint - position) * (coneHeight / 2));
	ofDrawSphere(spherePos, angle);
	// draw a cone object oriented towards aim position using the lookAt transformation
   // matrix.  The "up" vector is (0, 1, 0)
   //
	ofPushMatrix();
	glm::mat4 m = glm::lookAt(position, aimPoint, glm::vec3(0, 1, 0));
	ofMultMatrix(glm::inverse(m));

	ofRotateDeg(-90, 1, 0, 0);
	ofSetColor(ofColor::lightGray);
	ofDrawCone(angle, 1);
	ofPopMatrix();

	//draw aimpoint
	ofPushMatrix();
	ofSetColor(ofColor::purple);
	ofDrawSphere(aimPoint, aimtPointRadius);
	ofPopMatrix();

}


//--------------------------------------------------------------
void ofApp::setup() {
	//GUI setup
	gui.setup();
	gui.add(intensitySlider.setup("Light Intensity", 200, 1, 500));
	gui.add(powerSlider.setup("Phong Exponent", 500, 1, 10000));
	gui.add(ambientSlider.setup("Ambient Intensity", 0.05, 0, 0.5));
	gui.add(spotLightSlider.setup("Spot Light Intensity", 0.05, 0, 0.5));
	gui.add(shadowSlider.setup("Shadow Darkness", 0.05, 0, 0.5));
	//
	theCam = &mainCam;
	mainCam.setDistance(10);
	mainCam.setNearClip(.1);

	sideCam.setPosition(glm::vec3(10, 0, 0));
	sideCam.lookAt(glm::vec3(0, 0, 0));
	sideCam.setNearClip(.1);

	previewCam.setPosition(renderCam.position);
	previewCam.lookAt(renderCam.aim);
	previewCam.setNearClip(.1);
	//floor mesh
	Mesh *mesh = new Mesh();
	mesh->mVertices.clear();
	mesh->mVertices.push_back(Vertex(glm::vec3(-100, 10, -50)));//0
	mesh->mVertices.push_back(Vertex(glm::vec3(-100, -10, 50)));//1
	mesh->mVertices.push_back(Vertex(glm::vec3(100, 10, -50)));//2
	mesh->mVertices.push_back(Vertex(glm::vec3(100, -10, 50)));//3
	//mesh->mVertices.push_back(Vertex(glm::vec3(0, 4, 0)));//4
	//mesh->triList.clear();
	mesh->triList.push_back(Triangle({ 0, 2, 1 }));
	mesh->triList.push_back(Triangle({ 1, 2, 3 }));
	//mesh->triList.push_back(Triangle({ 1, 3, 4 }));
	//mesh->triList.push_back(Triangle({ 3, 2, 4 }));
	//mesh->triList.push_back(Triangle({ 2, 0, 4 }));
	//mesh->triList.push_back(Triangle({ 0, 1, 4 }));
	mesh->diffuseColor = ofColor::aquamarine;
	mesh->setReflective(0.9);
	scene.push_back(mesh);

	//textures
	//floorTexture.load("images/stone_texture_test.jpg");
	//floorSpecular.load("images/stone_floor_spec.jpg");
	skyTexture.load("images/skytexture.jpg");

	// each sphere has y = 1 has offset for their radius
	Sphere* s1 = new Sphere(glm::vec3(-6, 1, 10), 1, ofColor::blueSteel);
	s1->setReflective(0.5);
	scene.push_back(s1);
	Sphere* s2 = new Sphere(glm::vec3(4, 3, 7), 1, ofColor::aquamarine);
	s2->setReflective(0.5);
	scene.push_back(s2);
	Sphere* s3 = new Sphere(glm::vec3(0, 4, -5), 1, ofColor::aliceBlue);
	s3->setReflective(0.5);
	scene.push_back(s3);
	// Plane* floorPlane = new Plane(glm::vec3(3, 0, 2), glm::vec3(0, 1, 0), ofColor::skyBlue, 80, 80);
	/*floorPlane->plane.rotateDeg(10, 1, 0, 0);
	floorPlane->setTexture(floorTexture);
	floorPlane->setReflective(1);*/
	//floorPlane->setSpecular(floorSpecular);
	//scene.push_back(floorPlane);

	Plane* sky = new Plane(glm::vec3(3, 7, 2), glm::vec3(0, -1, 0), ofColor::lightSkyBlue, 200, 200);
	sky->noshading = true;
	sky->setTexture(skyTexture);
	scene.push_back(sky);

	//
	Light* mainLight = new Light(glm::vec3(-2, 2, 50), 1);
	lightList.push_back(mainLight);
	Light* subLight1 = new Light(glm::vec3(5, 0, 60), 1);
	lightList.push_back(subLight1);
	// Light* backLight1 = new Light(glm::vec3(-10, 20, -15), 1);
	// lightList.push_back(backLight1);
	//Light* backLight2 = new Light(glm::vec3(10, 20, -15), 1);
	//lightList.push_back(backLight2);
	//Light* backLight3 = new Light(glm::vec3(-10, 20, -200), 1);
	//lightList.push_back(backLight3);
	//Light* backLight4 = new Light(glm::vec3(10, 20, -200), 1);
	//lightList.push_back(backLight4);
	//SpotLight* spotLight1 = new SpotLight(glm::vec3(-10, 10, 3), glm::vec3(0, 0, 0), 0.2, 0.2, 1);
	//spotLights.push_back(spotLight1);
	//scene.push_back(spotLight1);
	//SpotLight* spotLight2 = new SpotLight(glm::vec3(4, 9, 7), glm::vec3(2, 0, 0), 0.2, 0.2, 1);
	//spotLights.push_back(spotLight2);
	//scene.push_back(spotLight2);
}

//--------------------------------------------------------------
void ofApp::update() {
	for (int i = 0; i < lightList.size(); i++) {
		lightList[i]->setIntensity(intensitySlider);
	}
}

//--------------------------------------------------------------
void ofApp::draw() {

	if (bShowImage) {
		ofTexture tt;
		tt.allocate(image.getPixels());
		tt.draw(0, 0);
	}
	else {
		theCam->begin();
		if (rayTest) {
			glm::vec3 p1 = glm::vec3(3, 1, 4);
			glm::vec3 p2 = glm::vec3(4, 2, 5);
			glm::vec3 p3 = glm::vec3(5, 1, 4);
			ofSetColor(ofColor::darkGoldenRod);
			ofDrawLine(p1, p2);
			ofDrawLine(p2, p3);
			ofDrawLine(p1, p3);

		}
		for (int i = 0; i < scene.size(); i++) {
			SceneObject* current = scene[i];
			ofSetColor(current->diffuseColor);
			current->draw();
		}
		for (int j = 0; j < lightList.size(); j++) {
			Light* l = lightList[j];
			ofSetColor(ofColor::black);
			l->draw();
		}
		theCam->end();
	}

	//ray trace
	if (bStartRayTrace) {
		ofSetColor(ofColor::white);
		std::thread t = std::thread(&ofApp::rayTrace, this);
		t.detach();
		bStartRayTrace = false;
		bShowImage = true;
	}

	gui.draw();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key) {
	switch (key) {
	case OF_KEY_F1:
		theCam = &mainCam;
		break;
	case OF_KEY_F2:
		theCam = &sideCam;
		break;
	case OF_KEY_F3:
		theCam = &previewCam;
		break;
	case 'u':
		renderCam.setPosition(mainCam.getPosition());
		break;
	case 't':
		rayTest = !rayTest;
		break;
	case 'r':
		bStartRayTrace = !bStartRayTrace;
		break;
	case 'R':
		bShowImage = !bShowImage;
		break;
	case 'p':
		bPhongShading = !bPhongShading;
		break;
	case 'c':
		if (mainCam.getMouseInputEnabled()) mainCam.disableMouseInput();
		else mainCam.enableMouseInput();
		break;
	case ',':
		for (int i = 0; i < spotLights.size(); i++) {
			spotLights[i]->angle = glm::min(10.0, spotLights[i]->angle + 0.1);
		}
		break;
	case '.':
		for (int j = 0; j < spotLights.size(); j++) {
			spotLights[j]->angle = glm::max(0.0, spotLights[j]->angle - 0.1);
		}
		break;
	default:
		break;
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key) {

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y) {

}


//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button) {
	if (mainCam.getMouseInputEnabled()) return;
	glm::vec3 origin = mainCam.getPosition();
	glm::vec3 mouseWorld = mainCam.screenToWorld(glm::vec3(x, y, 0));
	glm::vec3 mouseDir = glm::normalize(mouseWorld - origin);
	for (int i = 0; i < spotLights.size(); i++) {
		SpotLight* spotLight = spotLights[i];

		if (spotLights[i]->coneSelected) {

			glm::vec3 oldPos = spotLight->position;
			float distance;
			bool hit = glm::intersectRayPlane(origin, mouseDir, oldPos, mainCam.getZAxis(), distance);
			glm::vec3 mousePos = origin + distance * mouseDir;
			glm::vec3 delta = mousePos - mouseLastPos;

			oldPos += delta;
			spotLight->position = oldPos;

			mouseLastPos = mousePos;

		}
		else if (spotLights[i]->aimPointSelected) {
			glm::vec3 oldPos = spotLight->aimPoint;
			float distance;
			bool hit = glm::intersectRayPlane(origin, mouseDir, oldPos, mainCam.getZAxis(), distance);
			glm::vec3 mousePos = origin + distance * mouseDir;
			glm::vec3 delta = mousePos - mouseLastPos;

			oldPos += delta;
			spotLight->aimPoint = oldPos;
			mouseLastPos = mousePos;
		}
	}
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button) {
	if (mainCam.getMouseInputEnabled()) return;



	glm::vec3 origin = mainCam.getPosition();
	glm::vec3 mouseWorld = mainCam.screenToWorld(glm::vec3(x, y, 0));
	glm::vec3 mouseDir = glm::normalize(mouseWorld - origin);
	Ray ray = Ray(origin, mouseDir);
	glm::vec3 point, normal, targetPoint, targetNormal;

	if (rayTest) {
		glm::vec2 barycentric;
		float dis = 0;
		glm::vec3 p1 = glm::vec3(3, 1, 4);
		glm::vec3 p2 = glm::vec3(4, 2, 5);
		glm::vec3 p3 = glm::vec3(5, 1, 4);
		hitTriangle = glm::intersectRayTriangle(ray.p, ray.d, p1, p2, p3, barycentric, dis);
		if (hitTriangle) {

			glm::vec3 intersectPoint = ray.evalPoint(dis);
			printf("barycentric %f,%f\n ", barycentric.x, barycentric.y);
			printf("distance %f \n", dis);
			printf("hit test triangle at %f,%f,%f\n ", intersectPoint.x, intersectPoint.y, intersectPoint.z);

		}

	}

	//to be looped if more spotlights are added
	for (int i = 0; i < spotLights.size(); i++) {
		//SpotLight* spotLight = spotLights[i];
		bool hitSpotLight = glm::intersectRaySphere(ray.p, ray.d, spotLights[i]->position, spotLights[i]->angle, point, normal);
		bool hitAimPoint = glm::intersectRaySphere(ray.p, ray.d, spotLights[i]->aimPoint, spotLights[i]->aimtPointRadius, targetPoint, targetNormal);

		if (hitSpotLight) {
			float distance;
			bool hit = glm::intersectRayPlane(origin, mouseDir, spotLights[i]->position, mainCam.getZAxis(), distance);
			mouseDownPos = origin + distance * mouseDir;
			mouseLastPos = mouseDownPos;
			spotLights[i]->coneSelected = true;
		}
		else if (hitAimPoint) {
			float distance;
			bool hit = glm::intersectRayPlane(origin, mouseDir, spotLights[i]->aimPoint, mainCam.getZAxis(), distance);
			mouseDownPos = origin + distance * mouseDir;
			mouseLastPos = mouseDownPos;
			spotLights[i]->aimPointSelected = true;
		}
	}
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button) {
	for (int i = 0; i < spotLights.size(); i++) {
		spotLights[i]->coneSelected = false;
		spotLights[i]->aimPointSelected = false;
	}
}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y) {

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h) {

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg) {

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo) {
	for (int i = 0; i < dragInfo.files.size(); i++)
	{
		cout << dragInfo.files[i] << endl;//print dropped file path
	}
	objectFile.open(dragInfo.files[0]);
	Mesh* meshPtr = new Mesh();
	Mesh& mesh = *meshPtr;
	mesh.mVertices.clear();
	mesh.mVertices.push_back(Vertex(glm::vec3(0, 0, 0)));//filler vertex to match with .obj indexing starting from 1
	mesh.triList.clear();
	int vCter = 0;
	int fCter = 0;
	string currentLine;

	while (!objectFile.eof()) {
		getline(objectFile, currentLine);
		cout << "current line " + currentLine << endl;
		if (currentLine[0] == '#' || currentLine.size() < 1)
			continue;

		const char delim = ' ';
		vector<string> out;
		tokenize(currentLine, delim, out);

		//first char of first string in string array (eq: v 1.00 2.00 3.00)
		if (out[0] == "v" && out[0].size() == 1) {//checkif starts with v and is not vt nor vn
			vCter++;
			float x = stof(out[1]);
			float y = stof(out[2]);
			float z = stof(out[3]);
			Vertex v = Vertex(glm::vec3(x, y, z));//create vertex with x, y, z coordinates
			mesh.mVertices.push_back(v);
		}

		else if (out[0] == "f") {//eq. f 1/1/1 2/2/2 3/3/3
			fCter++;
			int xInd = stoi(out[1]);
			int yInd = stoi(out[2]);
			int zInd = stoi(out[3]);
			Triangle tri = Triangle({ xInd,yInd,zInd });//only get the first int, aka v (not vt or vn)

			mesh.triList.push_back(tri);
		}

	}
	objectFile.close();
	mesh.isLoaded = true;
	cout << "Vertices count: " + std::to_string(vCter) << endl;
	cout << "Faces count: " + std::to_string(fCter) << endl;
	//cout << "Mesh size: " + std::to_string(36 * vCter / 1000) + "kB" << endl;
	// Add mesh to scene
	mesh.diffuseColor = ofColor::saddleBrown;
	mesh.setReflective(0);
	scene.push_back(meshPtr);
	
}


//function from https://java2blog.com/split-string-space-cpp/
//to split string by delimitation char
//args: str to split, delimitation char, vector<string> as output
void ofApp::tokenize(string const& str, const char delim, vector<string>& out)
{
	// construct a stream from the string
	stringstream ss(str);

	string s;
	while (getline(ss, s, delim)) {
		out.push_back(s);
	}
}

ofColor ofApp::recursiveShade(const Ray &ray, int recursionDepth)
{
	float shortestDistance = std::numeric_limits<float>::max();
	glm::vec3 intersectPoint;
	glm::vec3 intersectNormal;
	int objIndex = -1;
	bool hit = false;
	
	for (int k = 0; k < scene.size(); k++) {
		SceneObject* s = scene[k];
		glm::vec3 point;
		glm::vec3 normal;
		
		if (s->intersect(ray, point, normal)) {
			// if obj is hit, check distance
			float distance = glm::distance(point, renderCam.view.position);
			if (distance < shortestDistance) {
				shortestDistance = distance;
				intersectPoint = point;
				intersectNormal = normal;
				hit = true;
				objIndex = k;
			}

		}

	}



	ofColor shading(0, 0, 0);

	if (hit) {
		if (scene[objIndex]->noshading)
			shading = scene[objIndex]->texture_lookup_wrap(intersectPoint, false);
		else {
			shading = phong(intersectPoint, intersectNormal, scene[objIndex]->specularColor, powerSlider, scene[objIndex]);
			if (recursionDepth > 0) {
				glm::vec3 rOutDirection = glm::normalize(ray.d) - (2 * intersectNormal * glm::dot(glm::normalize(ray.d), intersectNormal)); //refer to http://paulbourke.net/geometry/reflected/
				glm::vec3 rOutPos = intersectPoint + glm::normalize(rOutDirection) * 0.001;
				Ray rOut = Ray(rOutPos, rOutDirection);
				shading += scene[objIndex]->getReflective() * scene[objIndex]->specularColor * recursiveShade(rOut, recursionDepth - 1);
			}
		}
	}

	return shading;
}

//--------------------------------------------------------------
//ray tracing function
//
void ofApp::rayTrace() {

	//image.clear();
	image.allocate(imageWidth, imageHeight, OF_IMAGE_COLOR);
	image.setUseTexture(false);

	for (int j = image.getHeight(); j > 0; j--) {
		for (int i = 0; i < image.getWidth(); i++) {
			float u = (i + 0.5) / image.getWidth();
			float v = (j + 0.5) / image.getHeight();
			Ray ray = renderCam.getRay(u, v);
			ofColor shading = recursiveShade(ray, 2);
			image.setColor(i, image.getHeight() - j, shading);
		}
	}
	image.save("/images/output.png");
}

ofColor ofApp::ambient() {
	ofColor result = ofColor::black;
	result += ofColor::white * ambientSlider;
	return result;
}

ofColor ofApp::lambert(const glm::vec3& p, const glm::vec3& norm, SceneObject* obj) {
	ofColor result = ambient();

	for (int i = 0; i < lightList.size(); i++) {

		glm::vec3 p1 = p + (glm::normalize(norm) * 0.001);
		Ray shadRay = Ray(p1, glm::normalize(lightList[i]->position - p1));
		bool shad = false;

		for (int j = 0; j < scene.size(); j++) {
			SceneObject* current = scene[j];
			glm::vec3 shadPoint;
			glm::vec3 shadNormal;

			if (current->intersect(shadRay, shadPoint, shadNormal)) {
				shad = true;
				continue;
			}
		}
		if (shad)
			continue;
		else {
			glm::vec3 l = glm::normalize(lightList[i]->position - p);//vector to light source
			float r = glm::pow(glm::distance(p, lightList[i]->position), 2);
			ofColor diffuseFactor = obj->texture_lookup_wrap(p, false);
			result += diffuseFactor * (lightList[i]->intensity / r) * glm::max(float(0), glm::dot(norm, l));
		}

	}

	//account for spotlights
	for (int i = 0; i < spotLights.size(); i++) {
		SpotLight* spotLight = spotLights[i];
		Ray r = Ray(p, glm::normalize(spotLight->position - p));
		glm::vec3 point, normal;
		if (spotLight->intersect(r, point, normal))
			result += ofColor::white * spotLightSlider;
	}

	return result;
}

ofColor ofApp::phong(const glm::vec3& p, const glm::vec3& norm, const ofColor specular, float power, SceneObject* obj) {
	
	ofColor result = lambert(p, norm, obj);

	for (int i = 0; i < lightList.size(); i++) {

		glm::vec3 v = glm::normalize(renderCam.position - p);
		glm::vec3 l = glm::normalize(lightList[i]->position - p);//vector to light source
		glm::vec3 h = glm::normalize(l + v);
		float dotProduct = glm::dot(norm, h);
		float r = glm::pow(glm::distance(p, lightList[i]->position), 2);
		ofColor specularFactor = specular;
		if (obj->specularImage.isAllocated())
			specularFactor = obj->texture_lookup_wrap(p, true);
		result += specularFactor * (lightList[i]->intensity / r) * glm::max(float(0), glm::pow(dotProduct, power));

	}
	return result;
}


