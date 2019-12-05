#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <jni.h>

extern void square_cube_(
    int* i,
    int* isquare,
    int* icube
);

JNIEXPORT void JNICALL Java_com_example_stl_Stl_square_1cube
  (JNIEnv *env, jclass cl, jint i, jintArray isquare, jintArray icube)
{
    int i_int = i;
    jint *isquare_int = (*env)->GetIntArrayElements(env, isquare, 0);
    jint *icube_int = (*env)->GetIntArrayElements(env, icube, 0);
    square_cube_(&i_int, isquare_int, icube_int);

    (*env)->SetIntArrayRegion(env, isquare, 0, 1, isquare_int);
    (*env)->SetIntArrayRegion(env, icube, 0, 1, icube_int);

    (*env)->ReleaseIntArrayElements(env, isquare, isquare_int, JNI_ABORT);
    (*env)->ReleaseIntArrayElements(env, icube, icube_int, JNI_ABORT);
}

/*int main()
{
    int i = 2;
    printf("i = %d\n", i);
    int isquare = 0;
    int icube = 0;
    square_cube_(&i, &isquare, &icube);
    printf("isquare = %d\n", isquare);
    printf("icube = %d\n", icube);
    return 0;
}*/
