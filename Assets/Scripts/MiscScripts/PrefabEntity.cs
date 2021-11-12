using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Unity.Entities;
using Unity.Transforms;
using Unity.Rendering;

//[GenerateAuthoringComponent]
//public struct PrefabEntityComponent : IComponentData
//{
//    public Entity prefabEntity;
//}

public class PrefabEntity : MonoBehaviour, IDeclareReferencedPrefabs, IConvertGameObjectToEntity
{
    public static Entity prefabEntity;
    public GameObject prefabGameObject;

    public void Convert(Entity entity, EntityManager dstManager, GameObjectConversionSystem conversionSystem)
    {
        Entity prefabEntity = conversionSystem.GetPrimaryEntity(prefabGameObject);
        dstManager.AddComponents(prefabEntity, new ComponentTypes(new ComponentType[] { typeof(Parent), typeof(LocalToWorld) }));
        dstManager.RemoveComponent(prefabEntity, new ComponentTypes(new ComponentType[] { typeof(AmbientProbeTag), typeof(LinkedEntityGroup), typeof(PerInstanceCullingTag) }));
        PrefabEntity.prefabEntity = prefabEntity;
    }

    public void DeclareReferencedPrefabs(List<GameObject> referencedPrefabs)
    {
        referencedPrefabs.Add(prefabGameObject);
    }
}
