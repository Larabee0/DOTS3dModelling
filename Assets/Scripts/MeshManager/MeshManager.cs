using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Entities;
using Unity.Transforms;
using Unity.Collections;
using Unity.Rendering;
using Unity.Mathematics;
using System.Linq;
using Unity.Physics;

public class MeshManager : MonoBehaviour
{
    private EntityManager entityManager;
    private EntityArchetype MeshManagerArch;
    private Entity SpawnedMeshManagerArch;
    [SerializeField] private UnityEngine.Material[] materials;
    [SerializeField] private Mesh handleMesh;
    [SerializeField] private Mesh originalMesh;
    private Mesh clonedMesh;
    private bool isCloned = false;

    public float ScaleFactor = 0.025f;

    public static Mesh staticHandleMesh;
    public static UnityEngine.Material staticHandleMaterial;
    //[SerializeField] MeshData meshData;

    // Start is called before the first frame update
    private void Start()
    {
        staticHandleMesh = handleMesh;
        staticHandleMaterial = materials[1];
        float startTime = Time.realtimeSinceStartup;
        entityManager = World.DefaultGameObjectInjectionWorld.EntityManager;

        endSimulationEntityCommandBufferSystem = World.DefaultGameObjectInjectionWorld.GetOrCreateSystem<EndSimulationEntityCommandBufferSystem>();
        CloneMesh();
        CreateEntityArchetype();
        SpawnAndSetEntity();

        print("Pre Start up Time: " + (Time.realtimeSinceStartup - startTime) * 1000f + "ms");
    }
    //float startTime = Time.realtimeSinceStartup;
    //print("Start up Time: " + (Time.realtimeSinceStartup - startTime) * 1000f + "ms");

    private void CloneMesh()
    {
        if (!isCloned)
        {
            clonedMesh = new Mesh
            {
                name = name + "clone",
                vertices = originalMesh.vertices,
                triangles = originalMesh.triangles,
                normals = originalMesh.normals,
                uv = originalMesh.uv
            };
            isCloned = true;
        }
    }

    private void SpawnAndSetEntity()
    {
        SpawnedMeshManagerArch = entityManager.CreateEntity(MeshManagerArch);

        entityManager.SetComponentData(SpawnedMeshManagerArch, new Translation
        {
            Value = new float3(transform.position)
        });
        entityManager.SetComponentData(SpawnedMeshManagerArch, new Scale
        {
            Value = 1
        });
        entityManager.SetSharedComponentData(SpawnedMeshManagerArch, new RenderMesh
        {
            mesh = clonedMesh,
            material = materials[0]
        });

        entityManager.AddBuffer<Vertices>(SpawnedMeshManagerArch);
        entityManager.AddBuffer<Triangles>(SpawnedMeshManagerArch);
        entityManager.AddBuffer<SimplifiedVertices>(SpawnedMeshManagerArch);
        entityManager.AddBuffer<SimplifiedTriangles>(SpawnedMeshManagerArch);
        entityManager.AddBuffer<SimplifiedEdges>(SpawnedMeshManagerArch);
        float SetUpFirst = Time.realtimeSinceStartup;
        SetUpInitialData();
        entityManager.AddComponent(SpawnedMeshManagerArch, typeof(MeshManagerPreStart));
        print("Set Up Initial Data Rise Time: " + (Time.realtimeSinceStartup - SetUpFirst) * 1000f + "ms");
    }
    //private EntityManager entityManager;
    private EntityArchetype VertDataStoreArch;
    private EntityArchetype TriangleDataStoreArch;
    EntityArchetype EdgeDataStoreArch;
    private EndSimulationEntityCommandBufferSystem endSimulationEntityCommandBufferSystem;
    private EntityQueryDesc MeshManagerQuery;
    private EntityQueryDesc TriangleDataStoreQuery;
    private EntityQueryDesc EdgeDataStoreQuery;
    private EntityQueryDesc VertDataStoreQuery;
   
    private void NewStartUpSystem()
    {
        NewSystemPreStarter();
        NativeArray<int> VerticesToSimplifiedVerts;
        Dictionary<Vector3, List<VDRelatedVert>> relatedVertsDict;
        (relatedVertsDict, VerticesToSimplifiedVerts) = VertexDataPreperation(SpawnedMeshManagerArch);
        TriangleDataPreperation(SpawnedMeshManagerArch, relatedVertsDict, VerticesToSimplifiedVerts);
        VerticesToSimplifiedVerts.Dispose();
        EdgeDataPreperation(SpawnedMeshManagerArch);
        NativeArray<SimplifiedVertices> SimpVertArrayFinal = VertexFinalisingJob(SpawnedMeshManagerArch, relatedVertsDict);
        EdgeFinalisingJob(SpawnedMeshManagerArch, SimpVertArrayFinal);
        TriangleFinalisingJob(SpawnedMeshManagerArch, SimpVertArrayFinal);
    }

    private void NewSystemPreStarter()
    {
        VertDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(VertDataStore) }
        };
        EdgeDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(EdgeDataStore) }
        };
        TriangleDataStoreQuery = new EntityQueryDesc
        {
            All = new ComponentType[] { typeof(TriangleDataStore) }
        };
        MeshManagerQuery = new EntityQueryDesc { All = new ComponentType[] { typeof(MeshManagerType) } };
        VertDataStoreArch = entityManager.CreateArchetype(
            typeof(Translation),
            typeof(Scale),
            typeof(RenderMesh),
            typeof(LocalToWorld),
            typeof(LocalToParent),
            typeof(Parent),
            typeof(RenderBounds),
            typeof(PhysicsCollider),
            typeof(VertDataStore)
        );
        EdgeDataStoreArch = entityManager.CreateArchetype(
            typeof(Translation),
            typeof(NonUniformScale),
            typeof(Rotation),
            typeof(RenderMesh),
            typeof(LocalToWorld),
            typeof(LocalToParent),
            typeof(Parent),
            typeof(RenderBounds),
            typeof(PhysicsCollider),
            typeof(EdgeDataStore)
        );

        TriangleDataStoreArch = entityManager.CreateArchetype(
            typeof(Translation),
            typeof(Rotation),
            typeof(LocalToWorld),
            typeof(LocalToParent),
            typeof(Parent),
            typeof(PhysicsCollider),
            typeof(TriangleDataStore)
        );
    }
    private (Dictionary<Vector3, List<VDRelatedVert>>, NativeArray<int>) VertexDataPreperation(Entity MeshManagerReal)
    {
        DynamicBuffer<Vertices> VerticeBuffer = entityManager.GetBuffer<Vertices>(MeshManagerReal);
        DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
        NativeHashSet<Vector3> LocalVertsHashSet = new NativeHashSet<Vector3>(SimplifiedVerticesBuffer.Length, Allocator.Temp);
        foreach (Vertices vertex in VerticeBuffer)
        {
            LocalVertsHashSet.Add(vertex.LocalVertex);
        }
        NativeArray<Vector3> SimpVertTemp = LocalVertsHashSet.ToNativeArray(Allocator.Temp);
        LocalVertsHashSet.Dispose();

        int SimpVertTempLength = SimplifiedVerticesBuffer.Capacity = SimpVertTemp.Length;
        Dictionary<Vector3, List<VDRelatedVert>> relatedVertsDict = new Dictionary<Vector3, List<VDRelatedVert>>(SimpVertTempLength);
        for (int i = 0; i < SimpVertTemp.Length; i++)
        {
            SimplifiedVerticesBuffer.Add(new SimplifiedVertices { Vertex = SimpVertTemp[i] });
            relatedVertsDict.Add(SimpVertTemp[i], new List<VDRelatedVert>());
        }
        SimpVertTemp.Dispose();

        Dictionary<Vector3, PreAndPostStartUpSystem.SimpVertContainer> RelatedVertDict = new Dictionary<Vector3, PreAndPostStartUpSystem.SimpVertContainer>(SimpVertTempLength);
        for (int i = 0; i < SimplifiedVerticesBuffer.Length; i++)
        {
            RelatedVertDict.Add(SimplifiedVerticesBuffer[i].Vertex, new PreAndPostStartUpSystem.SimpVertContainer(new NativeList<int>(Allocator.Temp), i));
        }

        int VertBufferLength = VerticeBuffer.Length;
        NativeArray<int> VerticesToSimplifiedVerts = new NativeArray<int>(VertBufferLength, Allocator.Temp);
        Vector3 Temp;
        for (int i = 0; i < VertBufferLength; i++)
        {
            Temp = VerticeBuffer[i].LocalVertex;
            RelatedVertDict[Temp].RelatedVerts.Add(i);
            VerticesToSimplifiedVerts[i] = RelatedVertDict[Temp].Index;
            relatedVertsDict[Temp].Add(new VDRelatedVert { RelatedVert = i });
        }

        PreAndPostStartUpSystem.SimpVertContainer[] RelatedVertDictArray = RelatedVertDict.Values.ToArray();
        for (int i = 0; i < RelatedVertDictArray.Length; i++)
        {
            for (int k = 0; k < RelatedVertDictArray[i].RelatedVerts.Length; k++)
            {

                Vertices aVertex = VerticeBuffer[k];
                aVertex.SimplifiedVertex = RelatedVertDictArray[i].Index;
                VerticeBuffer[k] = aVertex;
            }
        }
        return (relatedVertsDict, VerticesToSimplifiedVerts);
    }

    private void TriangleDataPreperation(Entity MeshManagerReal, Dictionary<Vector3, List<VDRelatedVert>> relatedVertsDict, NativeArray<int> VerticesToSimplifiedVerts)
    {
        DynamicBuffer<Triangles> TrianglesBuffer = entityManager.GetBuffer<Triangles>(MeshManagerReal);
        int length = TrianglesBuffer.Length;
        DynamicBuffer<SimplifiedTriangles> SimplifiedTrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal);
        int Triangles = 0;
        for (int i = 0; i < length; i += 3)
        {

            Triangles aTriangle = TrianglesBuffer[i];
            aTriangle.SimplifiedTriangle = VerticesToSimplifiedVerts[aTriangle.GlobalTriangle];
            TrianglesBuffer[i] = aTriangle;

            Triangles bTriangle = TrianglesBuffer[i + 1];
            bTriangle.SimplifiedTriangle = VerticesToSimplifiedVerts[bTriangle.GlobalTriangle];
            TrianglesBuffer[i + 1] = bTriangle;


            Triangles cTriangle = TrianglesBuffer[i + 2];
            cTriangle.SimplifiedTriangle = VerticesToSimplifiedVerts[cTriangle.GlobalTriangle];
            TrianglesBuffer[i + 2] = cTriangle;
            SimplifiedTrianglesBuffer.Add(new SimplifiedTriangles
            {
                TriangleIndex = Triangles,
                SimpVert1Index = aTriangle.SimplifiedTriangle,
                SimpVert2Index = bTriangle.SimplifiedTriangle,
                SimpVert3Index = cTriangle.SimplifiedTriangle
            });
            Triangles++;
        }
    }

    private void EdgeDataPreperation(Entity MeshManagerReal)
    {
        DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
        DynamicBuffer<SimplifiedEdges> SimplifiedEdgesBuffer = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal);
        DynamicBuffer<Triangles> TrianglesBuffer = entityManager.GetBuffer<Triangles>(MeshManagerReal);
        NativeHashSet<Vector3> Edgepos = new NativeHashSet<Vector3>(SimplifiedVerticesBuffer.Length * 2, Allocator.Temp);
        for (int i = 0; i < TrianglesBuffer.Length; i += 3)
        {
            NativeArray<int> EdgeCollection = new NativeArray<int>(new int[] { TrianglesBuffer[i].SimplifiedTriangle, TrianglesBuffer[i + 1].SimplifiedTriangle }, Allocator.Temp);
            Vector3 TempV3 = (SimplifiedVerticesBuffer[EdgeCollection[0]].Vertex + SimplifiedVerticesBuffer[EdgeCollection[1]].Vertex) / 2;
            if (!Edgepos.Contains(TempV3))
            {
                Edgepos.Add(TempV3);
                SimplifiedEdgesBuffer.Add(new SimplifiedEdges { EdgeIndex = SimplifiedEdgesBuffer.Length, SimpVert1Index = EdgeCollection[0], SimpVert2Index = EdgeCollection[1], EdgePosition = TempV3 });
            }
            EdgeCollection.Dispose();
            EdgeCollection = new NativeArray<int>(new int[] { TrianglesBuffer[i + 1].SimplifiedTriangle, TrianglesBuffer[i + 2].SimplifiedTriangle }, Allocator.Temp);
            TempV3 = (SimplifiedVerticesBuffer[EdgeCollection[0]].Vertex + SimplifiedVerticesBuffer[EdgeCollection[1]].Vertex) / 2;
            if (!Edgepos.Contains(TempV3))
            {
                Edgepos.Add(TempV3);
                SimplifiedEdgesBuffer.Add(new SimplifiedEdges { EdgeIndex = SimplifiedEdgesBuffer.Length, SimpVert1Index = EdgeCollection[0], SimpVert2Index = EdgeCollection[1], EdgePosition = TempV3 });
            }
            EdgeCollection.Dispose();
            EdgeCollection = new NativeArray<int>(new int[] { TrianglesBuffer[i + 2].SimplifiedTriangle, TrianglesBuffer[i].SimplifiedTriangle }, Allocator.Temp);
            TempV3 = (SimplifiedVerticesBuffer[EdgeCollection[0]].Vertex + SimplifiedVerticesBuffer[EdgeCollection[1]].Vertex) / 2;
            if (!Edgepos.Contains(TempV3))
            {
                Edgepos.Add(TempV3);
                SimplifiedEdgesBuffer.Add(new SimplifiedEdges { EdgeIndex = SimplifiedEdgesBuffer.Length, SimpVert1Index = EdgeCollection[0], SimpVert2Index = EdgeCollection[1], EdgePosition = TempV3 });
            }

            EdgeCollection.Dispose();
        }
    }

    private NativeArray<SimplifiedVertices> VertexFinalisingJob(Entity MeshManagerReal, Dictionary<Vector3, List<VDRelatedVert>> relatedVertsDict)
    {
        DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
        int length = SimplifiedVerticesBuffer.Length;
        NativeArray<Entity> VertexEntities = entityManager.CreateEntity(VertDataStoreArch, length, Allocator.Temp);
        NativeHashMap<Entity, int> simpVertIndexMap = new NativeHashMap<Entity, int>(length, Allocator.TempJob);
        for (int i = 0; i < length; i++)
        {
            SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
            SimplifiedVertices simplifiedVertex = SimplifiedVerticesBuffer[i];
            simplifiedVertex.VertEntity = VertexEntities[i];
            SimplifiedVerticesBuffer[i] = simplifiedVertex;

            DynamicBuffer<Child> children = entityManager.GetBuffer<Child>(MeshManagerReal);
            children.Add(new Child { Value = VertexEntities[i] });
            simpVertIndexMap.Add(VertexEntities[i], i);


            SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal);
            NativeArray<VDRelatedVert> vDRelatedVerts = new NativeArray<VDRelatedVert>(relatedVertsDict[SimplifiedVerticesBuffer[i].Vertex].ToArray(), Allocator.Temp);
            DynamicBuffer<VDRelatedVert> RelatedVertsBuffer = entityManager.AddBuffer<VDRelatedVert>(VertexEntities[i]);
            RelatedVertsBuffer.CopyFrom(vDRelatedVerts);
            vDRelatedVerts.Dispose();
        }
        VertexEntities.Dispose();
        NativeArray<SimplifiedVertices> SimpVertArrayFinal = (SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal)).ToNativeArray(Allocator.Temp);
        //SpawnVertexes vertexJob = new SpawnVertexes
        //{
        //    VDSTypeHandle = entityManager.GetComponentTypeHandle<VertDataStore>(false),
        //    colliderTypeHandle = entityManager.GetComponentTypeHandle<PhysicsCollider>(false),
        //    translationTypeHandle = entityManager.GetComponentTypeHandle<Translation>(false),
        //    vertexTypeHandle = entityManager.GetEntityTypeHandle(),
        //    SimplifiedVertices = SimplifiedVerticesBuffer,
        //    entityCommandBuffer = endSimulationEntityCommandBufferSystem.CreateCommandBuffer().AsParallelWriter(),
        //    MeshManagerEntity = MeshManagerReal,
        //    simpVertIndexMap = simpVertIndexMap
        //};
        //
        //EntityQuery vertexQueryForJob = entityManager.CreateEntityQuery(VertDataStoreQuery);
        //simpVertIndexMap.Dispose(vertexJob.ScheduleParallel(vertexQueryForJob, 1)).Complete();
        return SimpVertArrayFinal;
    }

    private void EdgeFinalisingJob(Entity MeshManagerReal, NativeArray<SimplifiedVertices> SimpVertArrayFinal)
    {
        DynamicBuffer<SimplifiedEdges> SimplifiedEdgesBuffer = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal);
        int length = SimplifiedEdgesBuffer.Length;
        NativeArray<Entity> EdgeEntities = entityManager.CreateEntity(EdgeDataStoreArch, length, Allocator.Temp);
        NativeHashMap<Entity, int> simpEdgeIndexMap = new NativeHashMap<Entity, int>(length, Allocator.TempJob);
        for (int i = 0; i < length; i++)
        {
            SimplifiedEdgesBuffer = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal);
            SimplifiedEdges SimplifiedEdge = SimplifiedEdgesBuffer[i];
            SimplifiedEdge.EdgeEntity = EdgeEntities[i];
            SimplifiedEdge.SimpVert1Entity = SimpVertArrayFinal[SimplifiedEdge.SimpVert1Index].VertEntity;
            SimplifiedEdge.SimpVert2Entity = SimpVertArrayFinal[SimplifiedEdge.SimpVert2Index].VertEntity;
            SimplifiedEdgesBuffer[i] = SimplifiedEdge;

            DynamicBuffer<Child> children = entityManager.GetBuffer<Child>(MeshManagerReal);
            children.Add(new Child { Value = EdgeEntities[i] });
            simpEdgeIndexMap.Add(EdgeEntities[i], i);

        }
        EdgeEntities.Dispose();

        //SpawnEdges edgeJob = new SpawnEdges
        //{
        //    EDSTypeHandle = entityManager.GetComponentTypeHandle<EdgeDataStore>(false),
        //    colliderTypeHandle = entityManager.GetComponentTypeHandle<PhysicsCollider>(false),
        //    translationTypeHandle = entityManager.GetComponentTypeHandle<Translation>(false),
        //    scaleTypeHandle = entityManager.GetComponentTypeHandle<NonUniformScale>(false),
        //    RotationTypeHandle = entityManager.GetComponentTypeHandle<Rotation>(false),
        //    edgeTypeHandle = entityManager.GetEntityTypeHandle(),
        //    SimplifiedVertices = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
        //    entityCommandBuffer = endSimulationEntityCommandBufferSystem.CreateCommandBuffer().AsParallelWriter(),
        //    MeshManagerEntity = MeshManagerReal,
        //    SimplifiedEdges = entityManager.GetBuffer<SimplifiedEdges>(MeshManagerReal),
        //    simpEdgeIndexMap = simpEdgeIndexMap
        //};
        //
        //EntityQuery edgeQueryForJob = entityManager.CreateEntityQuery(EdgeDataStoreQuery);
        //simpEdgeIndexMap.Dispose(edgeJob.ScheduleParallel(edgeQueryForJob, 1)).Complete();
    }

    private void TriangleFinalisingJob(Entity MeshManagerReal,NativeArray<SimplifiedVertices> SimpVertArrayFinal)
    {
        DynamicBuffer<SimplifiedTriangles> SimplifiedTrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal);
        int length = SimplifiedTrianglesBuffer.Length;
        NativeArray<Entity> TriangleEntities = entityManager.CreateEntity(TriangleDataStoreArch, length, Allocator.Temp);
        NativeHashMap<Entity, int> simpTriangleIndexMap = new NativeHashMap<Entity, int>(length, Allocator.TempJob);
        for (int i = 0; i < length; i++)
        {
            SimplifiedTrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal);
            SimplifiedTriangles current = SimplifiedTrianglesBuffer[i];
            current.TriangleEntity = TriangleEntities[i];
            current.SimpVert1Entity = SimpVertArrayFinal[current.SimpVert1Index].VertEntity;
            current.SimpVert2Entity = SimpVertArrayFinal[current.SimpVert2Index].VertEntity;
            current.SimpVert3Entity = SimpVertArrayFinal[current.SimpVert3Index].VertEntity;

            DynamicBuffer<Child> children = entityManager.GetBuffer<Child>(MeshManagerReal);
            children.Add(new Child { Value = TriangleEntities[i] });
            simpTriangleIndexMap.Add(TriangleEntities[i], i);
        }
        TriangleEntities.Dispose();
        SimpVertArrayFinal.Dispose();

        //SpawnTriangles triangleJob = new SpawnTriangles
        //{
        //    TDSTypeHandle = entityManager.GetComponentTypeHandle<TriangleDataStore>(false),
        //    triangleTypeHandle = entityManager.GetEntityTypeHandle(),
        //    SimplifiedVertices = entityManager.GetBuffer<SimplifiedVertices>(MeshManagerReal),
        //    entityCommandBuffer = endSimulationEntityCommandBufferSystem.CreateCommandBuffer().AsParallelWriter(),
        //    MeshManagerEntity = MeshManagerReal,
        //    SimplifiedTriangles = entityManager.GetBuffer<SimplifiedTriangles>(MeshManagerReal),
        //    simpTriangleIndexMap = simpTriangleIndexMap
        //};
        //
        //EntityQuery triangleQueryForJob = entityManager.CreateEntityQuery(TriangleDataStoreQuery);
        //simpTriangleIndexMap.Dispose(triangleJob.ScheduleParallel(triangleQueryForJob, 1)).Complete();
    }

    private void OldStartUpSystem()
    {
        float SetUpFirst = Time.realtimeSinceStartup;
        OldSetUpVertexAndTriangleData();
        print("Set Up Vertex And Triangle Data Rise Time: " + (Time.realtimeSinceStartup - SetUpFirst) * 1000f + "ms");
        float SetUpThird = Time.realtimeSinceStartup;
        OldSetUpEdgeData();
        print("Set Up Edge Data Rise Time: " + (Time.realtimeSinceStartup - SetUpThird) * 1000f + "ms");
        print("Set Up Rise Time: " + (Time.realtimeSinceStartup - SetUpFirst) * 1000f + "ms");
        
        float SpawnFirst = Time.realtimeSinceStartup;
        OldSpawnVerticeIndicators();
        print("Spawn Vertice Indicators Rise Time: " + (Time.realtimeSinceStartup - SpawnFirst) * 1000f + "ms");
        float SpawnSecond = Time.realtimeSinceStartup;
        OldSpawnEdgeIndicators();
        print("Spawn Edge Indicators Rise Time: " + (Time.realtimeSinceStartup - SpawnSecond) * 1000f + "ms");
        float SpawnThird = Time.realtimeSinceStartup;
        OldSpawnTriangleIndicators();
        print("Spawn Triangle Indicators Rise Time: " + (Time.realtimeSinceStartup - SpawnThird) * 1000f + "ms");
        print("Spawn Rise Time: " + (Time.realtimeSinceStartup - SpawnFirst) * 1000f + "ms");
        
        float PostSpawn = Time.realtimeSinceStartup;
        PostSpawnDataSetting();
        print("Post Spawn Data Rise Time: " + (Time.realtimeSinceStartup - PostSpawn) * 1000f + "ms");
    }

    private void PostSpawnDataSetting()
    {
        DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(SpawnedMeshManagerArch);
        DynamicBuffer<SimplifiedEdges> SimplifiedEdgesBuffer = entityManager.GetBuffer<SimplifiedEdges>(SpawnedMeshManagerArch);
        DynamicBuffer<SimplifiedTriangles> SimplifiedTrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(SpawnedMeshManagerArch);
        NativeArray<SimplifiedVertices> SimplifiedVerticesArray = SimplifiedVerticesBuffer.ToNativeArray(Allocator.Temp);
        NativeArray<SimplifiedEdges> SimplifiedEdgesArray = SimplifiedEdgesBuffer.ToNativeArray(Allocator.Temp);
        NativeArray<SimplifiedTriangles> SimplifiedTrianglesArray = SimplifiedTrianglesBuffer.ToNativeArray(Allocator.Temp);
        //ConnectedEdges for vertices
        for (int i = 0; i < SimplifiedEdgesArray.Length; i++)
        {
            Entity Vert1Entitiy = SimplifiedVerticesBuffer[SimplifiedEdgesArray[i].SimpVert1Index].VertEntity;
            Entity Vert2Entitiy = SimplifiedVerticesBuffer[SimplifiedEdgesArray[i].SimpVert2Index].VertEntity;
            entityManager.GetBuffer<VDConnectedEdge>(Vert1Entitiy).Add(new VDConnectedEdge { ConnectedEdgeIndex = i, ConnectedEdgeEntity = SimplifiedEdgesArray[i].EdgeEntity });
            entityManager.GetBuffer<VDConnectedEdge>(Vert2Entitiy).Add(new VDConnectedEdge { ConnectedEdgeIndex = i, ConnectedEdgeEntity = SimplifiedEdgesArray[i].EdgeEntity });
            DynamicBuffer<VDConnectedVert> TempBuffer1 = entityManager.GetBuffer<VDConnectedVert>(Vert1Entitiy);
            DynamicBuffer<VDConnectedVert> TempBuffer2 = entityManager.GetBuffer<VDConnectedVert>(Vert2Entitiy);
            NativeArray<VDConnectedVert> VDConnectedVerts1Array = TempBuffer1.ToNativeArray(Allocator.Temp);
            NativeArray<VDConnectedVert> VDConnectedVerts2Array = TempBuffer2.ToNativeArray(Allocator.Temp);
            if (!VDConnectedVerts1Array.Contains(new VDConnectedVert { ConnectedVertIndex = SimplifiedEdgesArray[i].SimpVert2Index, ConnectedVertEntity = SimplifiedEdgesArray[i].SimpVert2Entity }))
            {
                TempBuffer1.Add(new VDConnectedVert { ConnectedVertIndex = SimplifiedEdgesArray[i].SimpVert2Index, ConnectedVertEntity = SimplifiedEdgesArray[i].SimpVert2Entity });
            }
            if (!VDConnectedVerts2Array.Contains(new VDConnectedVert { ConnectedVertIndex = SimplifiedEdgesArray[i].SimpVert1Index, ConnectedVertEntity = SimplifiedEdgesArray[i].SimpVert1Entity }))
            {
                TempBuffer2.Add(new VDConnectedVert { ConnectedVertIndex = SimplifiedEdgesArray[i].SimpVert1Index, ConnectedVertEntity = SimplifiedEdgesArray[i].SimpVert1Entity });
            }
            VDConnectedVerts1Array.Dispose();
            VDConnectedVerts2Array.Dispose();
        }
        //Connected vertices is already known.

        //ConnectedEdgesForTriangles
        for (int i = 0; i < SimplifiedTrianglesArray.Length; i++)
        {
            Entity Vert1Entitiy = SimplifiedVerticesBuffer[SimplifiedTrianglesArray[i].SimpVert1Index].VertEntity;
            Entity Vert2Entitiy = SimplifiedVerticesBuffer[SimplifiedTrianglesArray[i].SimpVert2Index].VertEntity;
            Entity Vert3Entitiy = SimplifiedVerticesBuffer[SimplifiedTrianglesArray[i].SimpVert3Index].VertEntity;

            entityManager.GetBuffer<VDConnectedTriangle>(Vert1Entitiy).Add(new VDConnectedTriangle { ConnectedTriangle = i, ConnectedTriangleEntity = SimplifiedTrianglesArray[i].TriangleEntity });
            entityManager.GetBuffer<VDConnectedTriangle>(Vert2Entitiy).Add(new VDConnectedTriangle { ConnectedTriangle = i, ConnectedTriangleEntity = SimplifiedTrianglesArray[i].TriangleEntity });
            entityManager.GetBuffer<VDConnectedTriangle>(Vert3Entitiy).Add(new VDConnectedTriangle { ConnectedTriangle = i, ConnectedTriangleEntity = SimplifiedTrianglesArray[i].TriangleEntity });
            
            DynamicBuffer<VDConnectedEdge> ConnectedEdgesVertexBuffer1 = entityManager.GetBuffer<VDConnectedEdge>(Vert1Entitiy);
            DynamicBuffer<VDConnectedEdge> ConnectedEdgesVertexBuffer2 = entityManager.GetBuffer<VDConnectedEdge>(Vert2Entitiy);
            DynamicBuffer<VDConnectedEdge> ConnectedEdgesVertexBuffer3 = entityManager.GetBuffer<VDConnectedEdge>(Vert3Entitiy);

            NativeArray<VDConnectedEdge> VDConnectedEdges1Array = ConnectedEdgesVertexBuffer1.ToNativeArray(Allocator.Temp);
            NativeArray<VDConnectedEdge> VDConnectedEdges2Array = ConnectedEdgesVertexBuffer2.ToNativeArray(Allocator.Temp);
            NativeArray<VDConnectedEdge> VDConnectedEdges3Array = ConnectedEdgesVertexBuffer3.ToNativeArray(Allocator.Temp);
            VDConnectedEdge[][] ConnectedEdges = new VDConnectedEdge[][] { VDConnectedEdges1Array.ToArray(), VDConnectedEdges2Array.ToArray(), VDConnectedEdges3Array.ToArray() };
            DynamicBuffer<TDConnectedEdges> TDConnectedEdges = entityManager.GetBuffer<TDConnectedEdges>(SimplifiedTrianglesArray[i].TriangleEntity);
            
            int secondaryK = 1;
            for (int k = 0; k < 3; k++)
            {
                if (secondaryK == 3)
                {
                    secondaryK = 0;
                }
                for (int h = 0; h < ConnectedEdges[k].Length; h++)
                {
                    for (int l = 0; l < ConnectedEdges[secondaryK].Length; l++)
                    {
                        if (ConnectedEdges[k][h].ConnectedEdgeIndex == ConnectedEdges[secondaryK][l].ConnectedEdgeIndex)
                        {
                            TDConnectedEdges.Add(new TDConnectedEdges { ConnectedEdgeIndex = ConnectedEdges[k][h].ConnectedEdgeIndex, ConnectedEdgeEntity = SimplifiedEdgesArray[ConnectedEdges[k][h].ConnectedEdgeIndex].EdgeEntity });
                            DynamicBuffer<EDConnectedTriangles> EDConnectedTrinalgesBuffer = entityManager.GetBuffer<EDConnectedTriangles>(SimplifiedEdgesArray[ConnectedEdges[k][h].ConnectedEdgeIndex].EdgeEntity);
                            EDConnectedTrinalgesBuffer.Add(new EDConnectedTriangles { ConnectedTriangleIndex = i, ConnectedTriangleEntity = SimplifiedTrianglesArray[i].TriangleEntity });
                            goto Next;
                        }
                    }
                }
                Next:
                secondaryK++;
            }
            VDConnectedEdges1Array.Dispose();
            VDConnectedEdges2Array.Dispose();
            VDConnectedEdges3Array.Dispose();
        }


        SimplifiedVerticesArray.Dispose();
        SimplifiedEdgesArray.Dispose();
        SimplifiedTrianglesArray.Dispose();
    }

    private void OldSpawnTriangleIndicators()
    {
        EntityArchetype TriangleDataStoreArch = entityManager.CreateArchetype(
               typeof(Translation),
               typeof(Rotation),
               typeof(LocalToWorld),
               typeof(LocalToParent),
               typeof(Parent),
               typeof(PhysicsCollider),
               typeof(TriangleDataStore)
        );

        DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(SpawnedMeshManagerArch);
        NativeArray<SimplifiedVertices> SimplifiedVerticesArray = SimplifiedVerticesBuffer.ToNativeArray(Allocator.Temp);
        DynamicBuffer<SimplifiedTriangles> TrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(SpawnedMeshManagerArch);
        int length = TrianglesBuffer.Length;
        NativeArray<Entity> TriangleEntities = entityManager.CreateEntity(TriangleDataStoreArch, length, Allocator.Temp);
        TrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(SpawnedMeshManagerArch);
        float SetUpSecond = Time.realtimeSinceStartup;
        for (int i = 0; i < length; i++)
        {
            Entity FreshSpawned = TriangleEntities[i];
            SimplifiedTriangles current = TrianglesBuffer[i];
            DynamicBuffer<Child> children = entityManager.GetBuffer<Child>(SpawnedMeshManagerArch);
            children.Add(new Child { Value = FreshSpawned });
            current.TriangleEntity = FreshSpawned;
            current.SimpVert1Entity = SimplifiedVerticesArray[current.SimpVert1Index].VertEntity;
            current.SimpVert2Entity = SimplifiedVerticesArray[current.SimpVert2Index].VertEntity;
            current.SimpVert3Entity = SimplifiedVerticesArray[current.SimpVert3Index].VertEntity;
            entityManager.SetComponentData(FreshSpawned, new Parent
            {
                Value = SpawnedMeshManagerArch
            });

            NativeArray<float3> vertices = new NativeArray<float3>(new float3[] { SimplifiedVerticesArray[current.SimpVert1Index].Vertex, SimplifiedVerticesArray[current.SimpVert2Index].Vertex, SimplifiedVerticesArray[current.SimpVert3Index].Vertex }, Allocator.Temp);
            NativeArray<int3> triangles = new NativeArray<int3>(new int3[] { new int3(0, 1, 2) }, Allocator.Temp);
            BlobAssetReference<Unity.Physics.Collider> colliderReference = Unity.Physics.MeshCollider.Create(vertices, triangles);
            triangles.Dispose();
            entityManager.SetComponentData(FreshSpawned, new PhysicsCollider
            {
                Value = colliderReference
            });
            //current.TrianglePosition = (vertices[0] + vertices[1] + vertices[2]) / 3;
            vertices.Dispose();
            entityManager.SetComponentData(FreshSpawned, new TriangleDataStore
            {
                TriangleIndex = i
            });


            entityManager.AddBuffer<TDConnectedEdges>(FreshSpawned);
            TrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(SpawnedMeshManagerArch);
            TrianglesBuffer[i] = current;
        }
        print("Triangle Entities Setting Time: " + (Time.realtimeSinceStartup - SetUpSecond) * 1000f + "ms");
        SimplifiedVerticesArray.Dispose();
        TriangleEntities.Dispose();

    }

    private void OldSpawnEdgeIndicators()
    {
        EntityArchetype EdgeDataStoreArch = entityManager.CreateArchetype(
               typeof(Translation),
               typeof(NonUniformScale),
               typeof(Rotation),
               typeof(RenderMesh),
               typeof(LocalToWorld),
               typeof(LocalToParent),
               typeof(Parent),
               typeof(UpdateVisualComponent),
               typeof(RenderBounds),
               typeof(PhysicsCollider),
               typeof(EdgeDataStore)
           );

        DynamicBuffer<SimplifiedEdges> SimplifiedEdges = entityManager.GetBuffer<SimplifiedEdges>(SpawnedMeshManagerArch);
        DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(SpawnedMeshManagerArch);
        NativeArray<SimplifiedVertices> SimplifiedVerticesArray = SimplifiedVerticesBuffer.ToNativeArray(Allocator.Temp);
        int length = SimplifiedEdges.Length;
        NativeArray<Entity> EdgeEntities = entityManager.CreateEntity(EdgeDataStoreArch, length, Allocator.Temp);        
        float SetUpFirst = Time.realtimeSinceStartup;
        for (int i = 0; i < length; i++)
        {
            Entity FreshSpawned = EdgeEntities[i];
            SimplifiedEdges current = entityManager.GetBuffer<SimplifiedEdges>(SpawnedMeshManagerArch)[i];
            DynamicBuffer<Child> children = entityManager.GetBuffer<Child>(SpawnedMeshManagerArch);
            children.Add(new Child { Value = FreshSpawned });
            current.EdgeEntity = FreshSpawned;
            current.SimpVert1Entity = SimplifiedVerticesArray[current.SimpVert1Index].VertEntity;
            current.SimpVert2Entity = SimplifiedVerticesArray[current.SimpVert2Index].VertEntity;

            entityManager.SetComponentData(FreshSpawned, new EdgeDataStore
            {
                EdgeIndex = i
            });
            entityManager.SetComponentData(FreshSpawned, new Parent
            {
                Value = SpawnedMeshManagerArch
            });
            
            entityManager.SetSharedComponentData(FreshSpawned, new RenderMesh
            {
                mesh = handleMesh,
                material = materials[1]
            });
            BlobAssetReference<Unity.Physics.Collider> collider = Unity.Physics.BoxCollider.Create(new BoxGeometry
            {
                Center = float3.zero,
                Orientation = quaternion.identity,
                Size = new float3(1, 1, 1)
            });
            entityManager.SetComponentData(FreshSpawned, new PhysicsCollider { Value = collider });
            entityManager.AddBuffer<EDConnectedTriangles>(FreshSpawned);
            SimplifiedEdges = entityManager.GetBuffer<SimplifiedEdges>(SpawnedMeshManagerArch);
            SimplifiedEdges[i] = current;
        }
        print("Edge Entities Setting Time: " + (Time.realtimeSinceStartup - SetUpFirst) * 1000f + "ms");
        SimplifiedVerticesArray.Dispose();
        EdgeEntities.Dispose();

    }

    private void OldSpawnVerticeIndicators()
    {
        EntityArchetype VertDataStoreArch = entityManager.CreateArchetype(
            typeof(Translation),
            typeof(Scale),
            typeof(RenderMesh),
            typeof(LocalToWorld),
            typeof(LocalToParent),
            typeof(Parent),
            typeof(RenderBounds),
            typeof(UpdateVisualComponent),
            typeof(PhysicsCollider),
            typeof(VertDataStore)
        );

        DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(SpawnedMeshManagerArch);
        int length = SimplifiedVerticesBuffer.Length;
        NativeArray<Entity> VertexEntities = entityManager.CreateEntity(VertDataStoreArch, length, Allocator.Temp);
        SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(SpawnedMeshManagerArch);
        float SetUpSecond = Time.realtimeSinceStartup;
        for (int i = 0; i < length; i++)
        {
            Entity FreshSpawned = VertexEntities[i];
            SimplifiedVertices current = SimplifiedVerticesBuffer[i];
            DynamicBuffer<Child> children = entityManager.GetBuffer<Child>(SpawnedMeshManagerArch);
            children.Add(new Child { Value = FreshSpawned });
            current.VertEntity = FreshSpawned;
            entityManager.SetSharedComponentData(FreshSpawned, new RenderMesh
            {
                mesh = handleMesh,
                material = materials[1]
            });
            entityManager.SetComponentData(FreshSpawned, new Parent
            {
                Value = SpawnedMeshManagerArch
            });
            entityManager.SetComponentData(FreshSpawned, new Scale
            {
                Value = 0.1f
            });
            BlobAssetReference<Unity.Physics.Collider> collider = Unity.Physics.BoxCollider.Create(new BoxGeometry
            {
                Center = float3.zero,
                Orientation = quaternion.identity,
                Size = new float3(1, 1, 1)
            });
            entityManager.SetComponentData(FreshSpawned, new PhysicsCollider { Value = collider });
            entityManager.SetComponentData(FreshSpawned, new VertDataStore
            {
                SimplifiedVertexIndex = i
            });
            entityManager.AddBuffer<VDConnectedVert>(FreshSpawned).Add(new VDConnectedVert { ConnectedVertIndex = i, ConnectedVertEntity = FreshSpawned });
            entityManager.AddBuffer<VDConnectedEdge>(FreshSpawned);
            entityManager.AddBuffer<VDConnectedTriangle>(FreshSpawned);
            NativeArray<VDRelatedVert> vDRelatedVerts = new NativeArray<VDRelatedVert>( relatedVertsDict[current.Vertex].ToArray(),Allocator.Temp);
            DynamicBuffer<VDRelatedVert> RelatedVertsBuffer = entityManager.AddBuffer<VDRelatedVert>(FreshSpawned);
            RelatedVertsBuffer.CopyFrom(vDRelatedVerts);
            
            vDRelatedVerts.Dispose();
            SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(SpawnedMeshManagerArch);
            SimplifiedVerticesBuffer[i] = current;
        }
        print("Vertex Entities Setting Time: " + (Time.realtimeSinceStartup - SetUpSecond) * 1000f + "ms");
        VertexEntities.Dispose();
        relatedVertsDict.Clear();
    }


    

    private readonly Dictionary<Vector3, List<VDRelatedVert>> relatedVertsDict = new Dictionary<Vector3, List<VDRelatedVert>>();
    private void OldSetUpVertexAndTriangleData()
    {
        DynamicBuffer<Vertices> VerticeBuffer = entityManager.GetBuffer<Vertices>(SpawnedMeshManagerArch);
        DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(SpawnedMeshManagerArch);
        HashSet<Vector3> LocalVertsHashSet = new HashSet<Vector3>();
        foreach (Vertices vertex in VerticeBuffer)
        {
            LocalVertsHashSet.Add(vertex.LocalVertex);
        }
        Vector3[] SimpVertTemp = LocalVertsHashSet.ToArray();
        SimplifiedVerticesBuffer.Capacity = SimpVertTemp.Length;
        for (int i = 0; i < SimpVertTemp.Length; i++)
        {
            SimplifiedVerticesBuffer.Add(new SimplifiedVertices { Vertex = SimpVertTemp[i] });
            relatedVertsDict.Add(SimpVertTemp[i], new List<VDRelatedVert>());
        }
        Dictionary<Vector3, SimpVertContainer> RelatedVertDict = new Dictionary<Vector3, SimpVertContainer>(SimplifiedVerticesBuffer.Capacity);
        for (int i = 0; i < SimplifiedVerticesBuffer.Length; i++)
        {
            RelatedVertDict.Add(SimplifiedVerticesBuffer[i].Vertex, new SimpVertContainer(new List<int>(), i));
        }
        int[] VerticesToSimplifiedVerts = new int[VerticeBuffer.Capacity];
        Vector3 Temp;
        for (int i = 0; i < VerticeBuffer.Length; i++)
        {
            Temp = VerticeBuffer[i].LocalVertex;
            RelatedVertDict[Temp].RelatedVerts.Add(i);
            VerticesToSimplifiedVerts[i] = RelatedVertDict[Temp].Index;
            relatedVertsDict[Temp].Add(new VDRelatedVert { RelatedVert = i });
        }
        RelatedVertDict.Values.ToList().ForEach(Value =>
        {
            Value.RelatedVerts.ForEach(RelatedVert =>
            {
                Vertices aVertex = VerticeBuffer[RelatedVert];
                aVertex.SimplifiedVertex = Value.Index;
                VerticeBuffer[RelatedVert] = aVertex;
            });
        });
        DynamicBuffer<Triangles> TrianglesBuffer = entityManager.GetBuffer<Triangles>(SpawnedMeshManagerArch);
        int length = TrianglesBuffer.Length;
        DynamicBuffer<SimplifiedTriangles> SimplifiedTrianglesBuffer = entityManager.GetBuffer<SimplifiedTriangles>(SpawnedMeshManagerArch);
        int Triangles = 0;
        for (int i = 0; i < length; i+=3)
        {

            Triangles aTriangle = TrianglesBuffer[i];
            aTriangle.SimplifiedTriangle = VerticesToSimplifiedVerts[aTriangle.GlobalTriangle];
            TrianglesBuffer[i] = aTriangle;

            Triangles bTriangle = TrianglesBuffer[i+1];
            bTriangle.SimplifiedTriangle = VerticesToSimplifiedVerts[bTriangle.GlobalTriangle];
            TrianglesBuffer[i + 1] = bTriangle;


            Triangles  cTriangle = TrianglesBuffer[i + 2];
            cTriangle.SimplifiedTriangle = VerticesToSimplifiedVerts[cTriangle.GlobalTriangle];
            TrianglesBuffer[i + 2] = cTriangle;
            SimplifiedTrianglesBuffer.Add(new SimplifiedTriangles
            {
                TriangleIndex = Triangles,
                SimpVert1Index = aTriangle.SimplifiedTriangle,
                SimpVert2Index = bTriangle.SimplifiedTriangle,
                SimpVert3Index = cTriangle.SimplifiedTriangle
            });
            Triangles++;
        }
    }
    private void OldSetUpEdgeData()
    {
        DynamicBuffer<Triangles> TrianglesBuffer = entityManager.GetBuffer<Triangles>(SpawnedMeshManagerArch);
        DynamicBuffer<SimplifiedVertices> SimplifiedVerticesBuffer = entityManager.GetBuffer<SimplifiedVertices>(SpawnedMeshManagerArch);
        DynamicBuffer<SimplifiedEdges> SimplifiedEdgesBuffer = entityManager.GetBuffer<SimplifiedEdges>(SpawnedMeshManagerArch);
        HashSet<Vector3> Edgepos = new HashSet<Vector3>();
        // This should be quick to run as we advance through SimplifiedTriangles at 3 indexes per loop.
        for (int i = 0; i < TrianglesBuffer.Length; i += 3)
        {
            int[] Temp = new int[] { TrianglesBuffer[i].SimplifiedTriangle, TrianglesBuffer[i + 1].SimplifiedTriangle };
            Vector3 TempV3 = (SimplifiedVerticesBuffer[Temp[0]].Vertex + SimplifiedVerticesBuffer[Temp[1]].Vertex) / 2;
            if (!Edgepos.Contains(TempV3))
            {
                Edgepos.Add(TempV3);
                SimplifiedEdgesBuffer.Add(new SimplifiedEdges { EdgeIndex = SimplifiedEdgesBuffer.Length, SimpVert1Index = Temp[0], SimpVert2Index = Temp[1], EdgePosition = TempV3 });
            }

            Temp = new int[] { TrianglesBuffer[i + 1].SimplifiedTriangle, TrianglesBuffer[i + 2].SimplifiedTriangle };
            TempV3 = (SimplifiedVerticesBuffer[Temp[0]].Vertex + SimplifiedVerticesBuffer[Temp[1]].Vertex) / 2;
            if (!Edgepos.Contains(TempV3))
            {
                Edgepos.Add(TempV3);
                SimplifiedEdgesBuffer.Add(new SimplifiedEdges { EdgeIndex = SimplifiedEdgesBuffer.Length, SimpVert1Index = Temp[0], SimpVert2Index = Temp[1], EdgePosition = TempV3 });
            }

            Temp = new int[] { TrianglesBuffer[i + 2].SimplifiedTriangle, TrianglesBuffer[i].SimplifiedTriangle };
            TempV3 = (SimplifiedVerticesBuffer[Temp[0]].Vertex + SimplifiedVerticesBuffer[Temp[1]].Vertex) / 2;
            if (!Edgepos.Contains(TempV3))
            {
                Edgepos.Add(TempV3);
                SimplifiedEdgesBuffer.Add(new SimplifiedEdges { EdgeIndex = SimplifiedEdgesBuffer.Length, SimpVert1Index = Temp[0], SimpVert2Index = Temp[1], EdgePosition = TempV3 });
            }
        }
    }

    public struct SimpVertContainer
    {
        public List<int> RelatedVerts { get; set; }
        public int Index { get; set; }
        public SimpVertContainer(List<int> relatedVerts, int index)
        {
            RelatedVerts = relatedVerts;
            Index = index;
        }
    }
    
    private void SetUpInitialData()
    {
        DynamicBuffer<Vertices> VerticeBuffer = entityManager.GetBuffer<Vertices>(SpawnedMeshManagerArch);
        DynamicBuffer<Triangles> TriangleBuffer = entityManager.GetBuffer<Triangles>(SpawnedMeshManagerArch);
        Vector3[] verts = clonedMesh.vertices;
        int length = VerticeBuffer.Capacity = verts.Length;
        for (int i = 0; i < length; i++)
        {
            VerticeBuffer.Add(new Vertices { GlobalVertex = transform.TransformPoint(verts[i]), LocalVertex = verts[i] });
        }
        int[] tris = clonedMesh.triangles;
        length = TriangleBuffer.Capacity = tris.Length;
        for (int i = 0; i < length; i++)
        {
            TriangleBuffer.Add(new Triangles { GlobalTriangle = tris[i] });
        }

    }

    private void CreateEntityArchetype()
    {
        MeshManagerArch = entityManager.CreateArchetype(
            typeof(Translation),
            typeof(Scale),
            typeof(RenderMesh),
            typeof(LocalToWorld),
            typeof(RenderBounds),
            typeof(MeshManagerType),
            typeof(Child)
        );
    }
    public float3[] ConvertLocalToGlobal(Vector3[] localVertices)
    {
        float3[] GlobalVertices = new float3[localVertices.Length];
        for (int i = 0; i < localVertices.Length; i++)
        {
            GlobalVertices[i]=localVertices[i];
        }
        return GlobalVertices;
    }
}

