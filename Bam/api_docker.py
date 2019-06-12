from os import path
import code
import docker as docker_sdk

HERE = path.abspath(path.dirname(__file__))
ROOT = path.dirname(HERE)
pth = path.join(HERE, 'Dockerfile')

image_name = "test05"
container_name = "test_05"
docker = docker_sdk.from_env()
image, logs = docker.images.build(dockerfile=pth, path=ROOT, tag=image_name)

created_image = image.tags[0]
if image_name in created_image:

    ports = [8080]
    port_bindings = {8080: 8080}

    docker.containers.run(image=image,
                          name=container_name,
                          volumes= {
                          'full_path_of_local_volume_to_bind': {
                          'bind': '/usr/local/src/app',
                          'mode': 'rw',
                          }
                          },
                          ports={8080: 8080},
                          detach=True)
    print('image create and container running')
else:
    print("Docker Didnot run")