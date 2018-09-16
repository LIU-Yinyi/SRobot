## SRobot Libraires
- Author: DWJ, Lin Xiaoya(All except SCity); Google(SCity Hash)
- Revised: Jackie Wu, blues, Rui Yang, hongwen000, Champion-Liu
- Date: 2018-9-16
- Version: 2.0.4
- Abstract: Pratical Libraries for Robot Development. Education Only.

---

### Outline
Basic class List:

```cpp
	class SScene;
	class SAerialRobot;
	class SRadian;
	class SDegree;
	class SAngle : public SRadian;
	class SEulerAngle;

	class SArena;
	class SGroundRobot;
	class SID;

	class SStream : public SString;
	class SStream;						//数据流
	class SRadian;						//弧度
	class SDegree;						//角度
	class SAngle;						//角
	class SQuaternion;					//四元数
	class SMatrix3;						//3*3矩阵
	class SMatrix4;						//4*4矩阵
	class SVector2;						//二维向量
	class SVector3;						//三维向量
	class SVector4;						//四维向量
	class SPoint4;						//四维向量
	class SLine2;						//平面线段
	class SLine3;						//空间线段
	class SLine4;						//四维空间线段
	class SRectangle2;					//平面矩形
	class STime;						//时间
	class SMath;						//数学库
	class SInputInterface;				//数据输入接口定义
	class SOutputInterface;				//数据输出接口定义
	class SMessage;						//消息
```
---

### Build and Install
For linux:

```bash
make -j16
sudo make install
```

---

### Use Git Submodule in Project
Add submodules in project like `Chickens`:

```bash
# create a common file directory
cd Chickens
mkdir submodules
cd submodules

# add submodules
git submodules add -b master --name submodules/SRobot https://github.com/Champion-Liu/SRobot.git
git commit -m "add submodules of SRobot"
git push
```


---

### Appendix: git submodules usages

##### Make changes, commit and checkout submodule files
Just go the submodule directory and use git as usual

##### List all currently configured submodules     
```
git submodule
```   
or   
```
git submodule status
```   

##### Show information about a submodule
```
git remote show <remote>
```

##### Add a new submodule
Beware of the submodule name you choose: If you use a forward slash (/) git will think you want to delete the submodule and want to add all the files in the submodule directory. Please DON'T use a forward slash after the submodule name.

1. Run `git submodule add -b <branch> --name <name> <repository-path-or-url>`
2. Add the `.gitmodule` file and submodule folder to the superproject index
3. Commit both files on the superproject

##### Remove a submodule
1. Delete the relevant line from the `.gitmodules` file
2. Delete the relevant section from `.git/config`
3. Run `git rm --cached <submodule-path>` (no trailing slash)
4. Commit the superproject
5. Delete the now untracked submodule files

##### Clone a project with submodules
1. Clone the superproject as usual
2. Run `git submodule init` to init the submodules
3. Run `git submodule update` to have the submodules on a detached HEAD   

or

Run `git clone --recurse-submodules ssh://user@domain.tld/repo.git`

##### See all changes on submodules
```
git diff --submodule
```

##### Update the submodules to the lastest changes on their respective branches
```
git submodule update --remote
```